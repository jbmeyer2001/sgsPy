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
	DALTableMetadata& tableMetadata,
	std::vector<RasterBandMetaData>& bands,
	int width,
	int height,
	int nComp)
{	

	//read full bands into p_data
	for (size_t i = 0; i < bands.size(); i++) {
		RasterBandIO(bands[i], bands[i].p_buffer, width, height, 1, 1, width, height, true)
	}										    //read = true

	//remove nodata values from p_data, ensuring data pixels are consecutive
	std::vector<T> noDataVals(bands.size());
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands.nan);
	}

	int nFeatures = 0;
	for (int i = 0; i < height * width; i++) {
		bool isNan = false;
		for (int b = 0; b < bandCount; b++) { 
			T val = reinterpret_cast<T *>(bands[b].p_buffer)[(y * width) + x]
			isNan = val == noDataVals[b];
			if (noData) {
				break;
			}
			reinterpret_cast<T *>(bands[b].p_buffer)[nFeatures] = val;
		}
		nFeatures += !isNan;
	}

	//create dal heterogen table containing data
	DALHeterogenTable table = DALHeterogenTable::empty(tableMetadata);
	DALTableFill(table, bands, nFeatures);

	//calculate pca
	const auto desc = oneapi::dal::pca::descriptor<float, Method>().set_component_count(nComp).set_deterministic(true);
       	const auto result = oneapi::dal::pca::train(desc);
	std::cout << "Eigenvectors: " << std::endl << result.get_eigenvectors() << std::endl;
	std::cout << "Eigenvalues: " << std::endl << result.get_eigenvalues() << std::endl;
	return;	
}

template <typename T>
void//dal::pca::train_result
calculatePCA(
	DALTableMetadata& tableMetadata,
	std::vector<RasterBandMetaData>& bands,
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

	const auto desc = oneapi::dal::pca::descriptor<float, Method>().set_component_count(nComp).set_deterministic(true);
	oneapi::dal::pca::partial_train_result<> partial_result;

	for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
		for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
			int xValid, yValid;
			bands[0].GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
		
			//read bands into memory	
			for (size_t i = 0; i < bands.size(); i++) {
				rasterBandIO(bands[i], bands[i].p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true);
			}													//read = true

			//remove nodata values
			int nFeatures;
			for (int x = 0; x < xValid; x++) {
				for (int y = 0; y < yValid; y++) {
					bool isNan = false;
					for (int b = 0; b < bandCount; b++) { 
						T val = reinterpret_cast<T *>(bands[i].p_buffer)[(y * xBlockSize) + x]
						isNan = val == noDataVals[b];
						if (noData) {
							break;
						}
						reinterpret_cast<T *>(bands[i].p_buffer)[nFeatures] = val;
					}
					nFeatures += !isNan;
				}
			}

			//calculate partial result
			DALHeterogenTable table = DALHeterogenTable::empty(tableMetadata);
			DALTableFill(table, bands, nFeatures);
			partial_result = oneapi::dal::partial_train(desc, partial_result, table);
		}
	}

	auto result = oneapi::dal::finalize_train(desc, partial_result);
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

	std::vector<dal_feature_type> featureTypes(bandCount, dal_numeric);
	std::vector<dal_data_type> dataTypes(bandCount);

	int xBlockSize, yBlockSize;
	p_raster->getRasterBand(0)->GetBlockSize(&xBlockSize, &yBlockSize);

	for (int i = 0; i < p_raster->getBandCount(); i++) {
		bands[i].p_band = p_raster->getRasterBand(i);
		bands[i].nan = bands[i].p_band->GetNoDataValue();
		bands[i].type = p_raster->getRasterBandType(i);
		bands[i].size = p_raster->getRasterBandTypeSize(i);
		bands[i].p_buffer = !isMEMDataset ?
			VSIMalloc3(xBlockSize, yBlockSize, bands[i].size) :
			VSIMalloc3(width, height, bands[i].size);
	
		dataTypes[i] = GDAL2DAL(bands[i].type);		
	}

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

	//create table metadata from (numeric) feature type and data types
	DALArray<DALFeatureType> featureTypesArr = DALArray<DALFeatureType>::wrap(featureTypes.data(), featureTypes.size());
	DALArray<DALDataType> dataTypesArr = DALArray<DALDataType>::wrap(dataTypes.data(), dataTypes.size());
	DALTableMetadata tableMetadata = DALTableMetadata(dataTypesArr, featureTypesArr);

	//get block size
	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	if (largeRaster) {
		calculatePCA(tableMetadata, bands, xBlockSize, yBlockSize, xBlocks, yBlocks, nComp);	
	}	
	else {
		calculatePCA(tableMetadata, bands, width, height, nComp);
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
