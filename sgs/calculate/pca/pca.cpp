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

template <typename T>
void //dal::pca::train_result 
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int bandCount,
	void *p_data,
	int width,
	int height,
	int nComp)
{
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
			size * static_cast<size_t>(bandCount),
			size * static_cast<size_t>(bandCount) * width
		);
	}

	//remove nodata values from p_data, ensuring data pixels are consecutive
	int index = 0;
	std::vector<T> noDataVals(bands.size());
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands.nan);
	}

	for (int i = 0; i < height * width; i++) {
		bool isNan = false;
		for (int b = 0; b < bandCount; b++) { 
			T val = reinterpret_cast<T *>(p_data)[((y * xBlockSize) + x) * bandCount + b]
			isNan = val == noDataVals[b];
			if (noData) {
				break;
			}
			reinterpret_cast<T *>(p_data)[index * bandCount + b] = val;
		}
		index += !isNan;
	}

	//create dal::homogen_table
	auto table = dal::homogen_table::wrap<T>(
		reinterpret_cast<T *>(p_data),
		index,
		bandCount,
		data_layout::row_major
	);

	//calculate pca
	const auto desc = dal::pca::descriptor<float, Method>().set_component_count(nComp).set_deterministic(true);
       	const auto result = dal::pca::train(desc);
	std::cout << "Eigenvectors: " << std::endl << result.get_eigenvectors() << std::endl;
	std::cout << "Eigenvalues: " << std::endl << result.get_eigenvalues() << std::endl;
	return;	
}

template <typename T>
void//dal::pca::train_result
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int bandCount,
	void *p_data,
	int xBlockSize,
	int yBlockSize,
	int xBlocks,
	int yBlocks
	int nComp)
{
	std::vector<T> noDataValues(bands.size());
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands.nan);
	}

	const auto desc = dal::pca::descriptor<float, Method>().set_component_count(nComp).set_deterministic(true);
	dal::pca::partial_train_result<> partial_result;

	for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
		for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
			int xValid, yValid;
			bands[0].GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
		
			//read bands into memory	
			for (size_t i = 0; i < bands.size(); i++) {
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
					size * static_cast<size_t>(bandCount),
					size * static_cast<size_t>(bandCount) * static_cast<size_t>(xBlockSize)		
				);
			}

			//remove nodata values
			int index;
			for (int x = 0; x < xValid; x++) {
				for (int y = 0; y < yValid; y++) {
					bool isNan = false;
					for (int b = 0; b < bandCount; b++) { 
						T val = reinterpret_cast<T *>(p_data)[((y * xBlockSize) + x) * bandCount + b]
						isNan = val == noDataVals[b];
						if (noData) {
							break;
						}
						reinterpret_cast<T *>(p_data)[index * bandCount + b] = val;
					}
					index += !isNan;
				}
			}

			//create dal::homogen_table and update partial result
			auto table = dal::homogen_table::wrap<T>(
				reinterpret_cast<T *>(p_data),
				index,
				bandCount,
				data_layout::row_major
			);

			//calculate partial result
			partial_result = dal::partial_train(desc, partial_result, table);
		}
	}

	auto result = dal::finalize_train(desc, partial_result);
	std::cout << "Eigenvectors: " << std::endl << result.get_eigenvectors() << std::endl;
	std::cout << "Eigenvalues: " << std::endl << result.get_eigenvalues() << std::endl;
	return;
}

GDALRasterWrapper *pca(
	GDALRasterWrapper *p_raster,
	size_t nComp,
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

	GDALDataType type = GDT_Int32;
	size_t size = sizeof(int);
	for (int i = 0; i < p_raster->getBandCount(); i++) {
		bands[i].p_band = p_raster->getRasterBand(i);
		bands[i].nan = bands[i].p_band->GetNoDataValue();

		GDALDataType bandType = p_raster->getRasterBandType(i);

		//ensure the most precise data type is the one used for all bands
		if (type != GDT_Float64 && bandType == GDT_Float64) {
			type = GDT_Float64;
			size = sizeof(double);
		}
		else if (type == GDT_Int32 && bandType == GDT_Float32) {
			type = GDT_Float32;
			size = sizeof(float);
		}
	}

	//allocate memory which will hold all bands
	int xBlockSize, yBlockSize;
	bands[0].p_band->GetBlockSize(&xBlockSize, &yBlockSize);
	
	if (isMEMDataset) {
		p_dataset = createVirtualDataset("MEM", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? double : float;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			addBandToMEMDataset(p_dataset, pcaBands[i]);
		}
	}
	else if (isVRTDataset){
		p_dataset = createVirtualDataset("VRT", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = float;
			pcaBands[i].size = sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			createVRTBandDataset(p_dataset, pcaBands[i], tempFolder, pcaBands[i].name, VRTBandInfo, driverOptions);
		}
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();

		if (extension == ".tif") {
			driver = "GTiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now.");
		}
		
		bool useTiles = xBlockSize != width &&
				yBlockSize != height;

		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = float;
			pcaBands[i].size = sizeof(float);
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

	void *p_data = isMEMDataset ?
		VSIMalloc3(height * width, numBlocks, size) :
		VSIMalloc3(xBlockSize * yBlockSize, numBlocks, size);	

	//for each block...
	
	//get block size of first band
	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	if (largeRaster) {
		calculatePCA(pcaBands, type, size, bandCount, p_data, xBlockSize, yBlockSize, xBlocks, yBlocks, nComp);	
	}	
	else {
		calculatePCA(pcaBands, type, size, bandCount , p_data, width, height, nComp);
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
