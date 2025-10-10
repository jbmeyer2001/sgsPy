/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of PCA
 * Author: Joseph Meyer
 * Date: October, 2025
 *
 ******************************************************************************/

#include <iostream>

#include "helper.h"
#include "raster.h"

/*
#include "oneapi/dal/algo/pca/detail/train_ops.hpp"
#include "oneapi/dal/partial_train.hpp"
#include "oneapi/dal/train.hpp"
#include "oneapi/dal/algo/pca/common.hpp"
*/

//#include "oneapi/dal.hpp"
#include "oneapi/dal.hpp"

typedef oneapi::dal::homogen_table	DALHomogenTable;

struct PCAResult {
	std::vector<std::vector<double>> eigenvectors;
	std::vector<double> eigenvalues;
};

template <typename T>
PCAResult
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int width,
	int height,
	int nComp)
{	
	int bandCount = static_cast<int>(bands.size());
	T *p_data = reinterpret_cast<T *>(VSIMalloc3(width * height, bandCount, size));

	std::vector<T> noDataVals(bandCount);
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}

	//read full bands into p_data
	for (size_t i = 0; i < bands.size(); i++) {
		bands[i].p_band->RasterIO(
			GF_Read,
			0,
			0,
			width,
			height,
			(void *)((size_t)p_data + i * size),	
			width,
			height,
			type,
			size * bands.size(),
			size * bands.size() * width
		);	
	}

	//remove nodata values from p_data, ensuring data pixels are consecutive
	int nFeatures = 0;
	for (int i = 0; i < height * width; i++) {
		bool isNan = false;
		for (int b = 0; b < bandCount; b++) { 
			T val = p_data[i * bandCount + b];
			isNan = val == noDataVals[b] || std::isnan(val);
			if (isNan) {
				break;
			}
			p_data[nFeatures * bandCount + b] = val;
		}
		nFeatures += !isNan;
	}
	
	//calculate pca
	DALHomogenTable table = DALHomogenTable::wrap<T>(p_data, nFeatures, bandCount, oneapi::dal::data_layout::row_major);
	const auto desc = oneapi::dal::pca::descriptor<float, oneapi::dal::pca::method::cov>().set_component_count(nComp).set_deterministic(true);
       	const auto result = oneapi::dal::train(desc, table);

	VSIFree(p_data);

	PCAResult retval;
	auto eigenvectors = result.get_eigenvectors();
	auto eigenvalues = result.get_eigenvalues();
	int64_t eigRows = eigenvectors.get_row_count();
	int64_t eigCols = eigenvectors.get_column_count();

	oneapi::dal::row_accessor<const float> eigVecAcc {eigenvectors};
	auto eigVecBlock = eigVecAcc.pull({0, eigRows});

	oneapi::dal::row_accessor<const float> eigValAcc {eigenvalues};
	auto eigValBlock = eigValAcc.pull({0, 1});

	retval.eigenvectors.resize(eigRows);
	retval.eigenvalues.resize(eigRows);
	for (int64_t i = 0; i < eigRows; i++) {
		retval.eigenvalues[i] = static_cast<double>(eigValBlock[i]);

		retval.eigenvectors[i].resize(eigCols);
		for (int64_t j = 0; j < eigCols; j++) {
			retval.eigenvectors[i][j] = static_cast<double>(eigVecBlock[i * eigCols + j]);
		}
	}
	
	return retval;	
}

template <typename T>
PCAResult
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int xBlockSize,
	int yBlockSize,
	int xBlocks,
	int yBlocks,
	int nComp)
{
	int bandCount = static_cast<int>(bands.size());
	T *p_data = reinterpret_cast<T *>(VSIMalloc3(xBlockSize * yBlockSize, bandCount, size));

	std::vector<T> noDataVals(bandCount);
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}

	const auto desc = oneapi::dal::pca::descriptor<float, oneapi::dal::pca::method::cov>().set_component_count(nComp).set_deterministic(true);
	oneapi::dal::pca::partial_train_result<> partial_result;

	for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
		for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
			int xValid, yValid;
			bands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
		
			//read bands into memory	
			for (size_t i = 0; i < bandCount; i++) {
				bands[i].p_band->RasterIO(
					GF_Read,
					xBlock * xBlockSize,
					yBlock * yBlockSize,
					xValid,
					yValid,
					(void *)((size_t)p_data + i * size),
					xValid,
					yValid,
					type,
					size * bands.size(),
					size * bands.size() * static_cast<size_t>(xBlockSize)
				);				
			}

			//remove nodata values
			int nFeatures;
			for (int x = 0; x < xValid; x++) {
				for (int y = 0; y < yValid; y++) {
					bool isNan = false;
					for (int b = 0; b < bandCount; b++) { 
						T val = p_data[((y * xBlockSize) + x) * bandCount + b];
						isNan = std::isnan(val) || val == noDataVals[b];
						if (isNan) {
							break;
						}
						p_data[nFeatures * bandCount + b] = val;
					}
					nFeatures += !isNan;
				}
			}

			//calculate partial result
			DALHomogenTable table = DALHomogenTable(p_data, nFeatures, bandCount, [](const T*){}, oneapi::dal::data_layout::row_major);
			partial_result = oneapi::dal::partial_train(desc, partial_result, table);
		}
	}

	auto result = oneapi::dal::finalize_train(desc, partial_result);
	
	VSIFree(p_data);

	PCAResult retval;
	auto eigenvectors = result.get_eigenvectors();
	auto eigenvalues = result.get_eigenvalues();
	int64_t eigRows = eigenvectors.get_row_count();
	int64_t eigCols = eigenvectors.get_column_count();

	oneapi::dal::row_accessor<const float> eigVecAcc {eigenvectors};
	auto eigVecBlock = eigVecAcc.pull({0, eigRows});

	oneapi::dal::row_accessor<const float> eigValAcc {eigenvalues};
	auto eigValBlock = eigValAcc.pull({0, 1});

	retval.eigenvectors.resize(eigRows);
	retval.eigenvalues.resize(eigRows);
	for (int64_t i = 0; i < eigRows; i++) {
		retval.eigenvalues[i] = static_cast<double>(eigValBlock[i]);

		retval.eigenvectors[i].resize(eigCols);
		for (int64_t j = 0; j < eigCols; j++) {
			retval.eigenvectors[i][j] = static_cast<double>(eigValBlock[i * eigCols + j]);
		}
	}
	
	return retval;	
}

GDALRasterWrapper *pca(
	GDALRasterWrapper *p_raster,
	int nComp,
	bool plot,
	bool largeRaster,
	std::string tempFolder,
	std::string filename,
	std::map<std::string, std::string> driverOptions)
{
	GDALAllRegister();

	int bandCount = p_raster->getBandCount();
	int height = p_raster->getHeight();
	int width = p_raster->getWidth();
	double *geotransform = p_raster->getGeotransform();
	std::string projection = std::string(p_raster->getDataset()->GetProjectionRef());
	
	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";
	GDALDataset *p_dataset = nullptr;
	
	std::vector<RasterBandMetaData> bands(bandCount);
	std::vector<RasterBandMetaData> pcaBands(nComp);
	std::vector<VRTBandDatasetInfo> VRTBandInfo;

	int xBlockSize, yBlockSize;
	p_raster->getRasterBand(0)->GetBlockSize(&xBlockSize, &yBlockSize);

	GDALDataType type = GDT_Float32;
	size_t size = sizeof(float);
	for (int i = 0; i < p_raster->getBandCount(); i++) {
		bands[i].p_band = p_raster->getRasterBand(i);
		bands[i].nan = bands[i].p_band->GetNoDataValue();

		if (p_raster->getRasterBandType(i) == GDT_Float64) {
			type = GDT_Float64;
			size = sizeof(double);
		}
	}

	if (isMEMDataset) {
		p_dataset = createVirtualDataset("MEM", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			addBandToMEMDataset(p_dataset, pcaBands[i]);
		}
	}
	else if (isVRTDataset){
		p_dataset = createVirtualDataset("VRT", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			createVRTBandDataset(p_dataset, pcaBands[i], tempFolder, pcaBands[i].name, VRTBandInfo, driverOptions);
		}
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();
		std::string driver;

		if (extension == ".tif") {
			driver = "GTiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now.");
		}
		
		bool useTiles = xBlockSize != width &&
				yBlockSize != height;

		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			if (useTiles) {
				pcaBands[i].xBlockSize = xBlockSize;
				pcaBands[i].yBlockSize = yBlockSize;
			}
		}

		p_dataset = createDataset(
			filename,
			driver,
			width,
			height, 
			geotransform,
			projection,
			pcaBands.data(),
			pcaBands.size(),
			useTiles,
			driverOptions
		);
	}

	//get block size
	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	PCAResult result;
	switch(type) {
		case GDT_Float32:
			result = largeRaster ?
				calculatePCA<float>(bands, type, size, xBlockSize, yBlockSize, xBlocks, yBlocks, nComp) :
				calculatePCA<float>(bands, type, size, width, height, nComp);
			break;
		case GDT_Float64:
			result = largeRaster ?
				calculatePCA<double>(bands, type, size, xBlockSize, yBlockSize, xBlocks, yBlocks, nComp) :
				calculatePCA<double>(bands, type, size, width, height, nComp);
			break;
		default:
			throw std::runtime_error("should not be here! GDALDataType should be one of Int32/Float32/Float64!");
	}

	std::cout << "got result." << std::endl;

	std::cout << "eigenvalues:" << std::endl;
	for (const double& val : result.eigenvalues) {
		std::cout << val << " ";
	}
	std::cout << std::endl;
	std::cout << "eigenvectors: " << std::endl;
	for (const std::vector<double>& vector : result.eigenvectors) {
		for (const double& val : vector) {
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}

	throw std::runtime_error("debugging stop.");
	return nullptr;
	//call rasterio to read as a specific data type in
	//	row or column major
	
	//iterate through checking for nan values, deleting those which
	//	are
	
	//create a dal::homogen_table with the data in-placde
	
	//incrementally calculate pca

	//end for each block...

	//get final pca calculation
	
	//for each block...
	
	//calculate output val of each principal component per-pixel
	
	//write values to output dataset
	
	//end for each block...
	
	//return
}

PYBIND11_MODULE(pca, m) {
	m.def("pca_cpp", &pca);
}
