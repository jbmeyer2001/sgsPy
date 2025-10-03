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
dal::pca::train_result 
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int bandCount,
	void *p_data,
	int height,
	int width,
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
		bool noData = false;
		for (int band = 0; band < bandCount; band++) {
			T val = getPixelValueDependingOnType<T>(type, p_data, i * bandCount + band);
			noData = val == noDataVals[i];
			if (noData) {
				break;
			}
			static_cast<T *>(p_data)[index * bandCount + band] = static_cast<T *>(p_data)[i * bandCount + band];
		}
		index += !noData;
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
       	return dal::train::(desc, table);	
}

template <typename T>
dal::pca::train_result
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int bandCount,
	void *p_data,
	int xBlockSize,
	int yBlockSize,
	int xBlocks,
	int yBlocks)
{

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

	//TODO Move this to new function
	dal::pca::partial_train_result<> partial_result;
	
	void *p_data = isMEMDataset ?
		VSIMalloc3(height * width, numBlocks, size) :
		VSIMalloc3(xBlockSize * yBlockSize, numBlocks, size);	

	//for each block...
	
	//get block size of first band
	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			bands[0].getActualBlockSize(xBlock, yBlock, xValid, yValid);
			
			//read blocks into memory
			for (size_t band = 0; band < bands.size(); band++) {
				band->RasterIO(
					GF_Read,
					xBlock * xBlockSize,
					yBlock * yBlockSize,
					xValid,
					yValid,
					(void *)((size_t)p_data + band * size),
					xValid,
					yValid,
					type,
					size * static_cast<size_t>(bandCount),
					size * static_cast<size_t>(bandCount) * xBlockSize,
					nullptr
				);
			
				int rows = 0; //number of data values (rows) where columns are number of features
							
				switch(type) {
					case GDT_Int32:
						rows = removeNoDataFromBlock<int>(noDataValues, p_data, type, bandCount, xValid, yValid);
						auto table = dal::homogen_table::wrap<int>(
							reinterpret_cast<int *>(p_data), 
							rows, 
							bandCount,
							data_layout::row_major);
						break;
					case GDT_Float32;
						rows = removeNoDataFromBlock<float>(noDataValues, p_data, type, bandCount, xValid, yValid);
						break;
					default: //can only be GDT_Flaot64
						rows = removeNoDataFromBlock<double>(noDataValues, p_data, type, bandCount, xValid, yValid);
				}
			}
		}
	}

	auto table = dal::homogen_table::wrap()
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
