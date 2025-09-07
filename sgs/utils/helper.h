/******************************************************************************
 *
 * Project: sgs
 * Purpose: helper functions
 * Author: Joseph Meyer
 * Date: August, 2025
 *
 ******************************************************************************/

#pragma once

#include <iostream>
#include <filesystem>
#include <gdal_priv.h>

#define MAXINT8		127
#define MAXINT16	32767

/**
 *
 */
struct RasterBandMetaData {
	GDALRasterBand *p_band = nullptr;
	void *p_buffer = nullptr;
	GDALDataType type = GDT_Unknown;
	size_t size = 0;
	std::string name = "";
	double nan = -1;
	int xBlockSize = -1;
	int yBlockSize = -1;
	std::mutex mutex;
};

/**
 *
 */
struct VRTBandDatasetInfo {
	GDALDataset *p_dataset = nullptr;
	std::string filename = "";

};

/**
 * Helper function for setting the pixel type of a particular
 * strat raster band, by determining the maximum strata value
 * and finding the smallest signed integer type which will
 * hold that value without overflow.
 *
 * The function also adds the type onto the end of a vector
 * of GDALDataType type so that the caller knows the type
 * for future reference.
 *
 * @param size_t maxStrata
 * @param std::vector<GDALDataType>& stratBandTypes
 * @returns size_t pixel type size
 */
inline size_t
setStratBandType(
	size_t maxStrata,
	std::vector<GDALDataType>& stratBandTypes
) {
	size_t pixelTypeSize;

	if (maxStrata <= MAXINT8) {
		stratBandTypes.push_back(GDT_Int8);
		pixelTypeSize = sizeof(int8_t);
	}
	else if (maxStrata <= MAXINT16) {
		stratBandTypes.push_back(GDT_Int16);
		pixelTypeSize = sizeof(int16_t);
	}
	else {
		stratBandTypes.push_back(GDT_Int32);
		pixelTypeSize = sizeof(int32_t);
	}

	return pixelTypeSize;
}

/**
 *
 */
inline void
setStratBandTypeAndSize(
	size_t maxStrata,
	GDALDataType *p_type,
	size_t *p_size
) {
	if (maxStrata <= MAXINT8) {
		*p_type = GDT_Int8;
		*p_size = sizeof(int8_t);
	}
	else if (maxStrata <= MAXINT16) {
		*p_type = GDT_Int16;
		*p_size = sizeof(int16_t);
	}
	else {
		*p_type = GDT_Int32;
		*p_size = sizeof(int32_t);
	}
}

/**
 * Helper function for reading a particular pixel from a raster
 * data buffer. The data buffer is cast to the type corresponding
 * to the GDALDataType type parameter, then indexed. The resulting
 * value is then cast to the double type.
 *
 * @param GDALDataType type 
 * @param void *p_data
 * @param size_t index
 * @returns double pixel val
 */
template <typename T>
inline T
getPixelValueDependingOnType(
	GDALDataType type,
	void *p_data,
	size_t index
) {
	switch (type) {
		case GDT_Int8:
			return static_cast<T>(((int8_t *)p_data)[index]);
		case GDT_UInt16:
			return static_cast<T>(((uint16_t *)p_data)[index]);
		case GDT_Int16:
			return static_cast<T>(((int16_t *)p_data)[index]);
		case GDT_UInt32:
			return static_cast<T>(((uint32_t *)p_data)[index]);
		case GDT_Int32:
			return static_cast<T>(((int32_t *)p_data)[index]);
		case GDT_Float32:
			return static_cast<T>(((float *)p_data)[index]);
		case GDT_Float64:
			return static_cast<T>(((double *)p_data)[index]);
		default:
			throw std::runtime_error("raster pixel data type not supported.");
	}
}

/**
 * Helper function for writing a particular pixel value to a strat
 * raster data buffer. The data buffer is cast to the type
 * corresponding to the GDALDataType type parameter. Then,
 * the value at the index provided is set to either 'strata' if isNan
 * is false, or -1 if isNan is true in the type required.
 *
 * @param GDALDataType type
 * @param void *p_data
 * size_t index
 * bool isNan
 * size_t strata
 */
inline void
setStrataPixelDependingOnType(
	GDALDataType type,
	void *p_data,
	size_t index,
	bool isNan,
	size_t strata
) {
	switch(type) {
		case GDT_Int8:
			reinterpret_cast<int8_t *>(p_data)[index] = isNan ?
				static_cast<int8_t>(-1) :
				static_cast<int8_t>(strata);
			break;
		case GDT_Int16: 
			reinterpret_cast<int16_t *>(p_data)[index] = isNan ?
				static_cast<int16_t>(-1) :
				static_cast<int16_t>(strata);
			break;
		case GDT_Int32: 
			reinterpret_cast<int32_t *>(p_data)[index] = isNan ?
				static_cast<int32_t>(-1) :
				static_cast<int32_t>(strata);
			break;
		default:
			throw std::runtime_error("strata pixel data type not supported.");
	}
}

/**
 * Helper function which prints a warning to the user if conversion
 * from the raster data type they're using for a strat raster
 * may result in errors in conversion to a 32 bit signed integer
 * type.
 *
 * @param GDALDataType type
 */
inline void
printTypeWarningsForInt32Conversion(GDALDataType type) {
	switch(type) {
		case GDT_UInt32:
			std::cout << "**warning** the pixel type of one of the bands given is an unsigned 32 bit integer. This may result in undefined behavior if the value is not castable to a 32-bit signed integer type." << std::endl;
			break;
		case GDT_Float32:
			std::cout << "**warning** the pixel type of one of the bands given is a 32 bit floating point value. This may result in undefined behavior if the value is not castable to a 32-bit signed integer type." << std::endl;
			break;
		case GDT_Float64:
			std::cout << "**warning** the pixel type of one of the bands given is a 64 bit floating point value. This may result in undefined behavior if the value is not castable to a 32-bit signed integer type." << std::endl;
		default:
			//don't care if converting the type to int32 won't result in overflow problems
			break;
	}
}

/**
 *
 */
inline GDALDataset *
createVirtualDataset(
	std::string driverName, 
	int width, 
	int height, 
	double *geotransform, 
	std::string projection) 
{
	GDALAllRegister();
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName(driverName.c_str());
	if (!p_driver) {
		throw std::runtime_error("unable to find dataset driver.");
	}

	GDALDataset *p_dataset = p_driver->Create(
		"",
		width,
		height,
		0,
		GDT_Unknown,
		nullptr
	);
	if (!p_dataset) {
		throw std::runtime_error("unable to create dataset with driver.");
	}

	CPLErr err = p_dataset->SetGeoTransform(geotransform);
	if (err) {
		throw std::runtime_error("error setting geotransform.");
	}

	err = p_dataset->SetProjection(projection.c_str());
	if (err) {
		throw std::runtime_error("error setting projection.");
	}

	return p_dataset;
}

/**
 *
 */
inline GDALDataset *
createDataset(
	std::string filename,
	std::string driverName, 
	int width, 
	int height, 
	double *geotransform, 
	std::string projection,
	RasterBandMetaData *bands,
	size_t bandCount,
	bool useTiles) 
{
	GDALAllRegister();
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName(driverName.c_str());
	if (!p_driver) {
		throw std::runtime_error("unable to find dataset driver.");
	}

	char **papszOptions = nullptr;
	if (useTiles) {
		const char *xBlockSizeOption = std::to_string(bands[0].xBlockSize).c_str();
		const char *yBlockSizeOption = std::to_string(bands[0].yBlockSize).c_str();
		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		papszOptions = CSLSetNameValue(papszOptions, "BLOCKXSIZE", xBlockSizeOption);
		papszOptions = CSLSetNameValue(papszOptions, "BLOCKYSIZE", yBlockSizeOption);
	}

	GDALDataset *p_dataset = p_driver->Create(
		filename.c_str(),
		width,
		height,
		bandCount,
		bands[0].type,
		papszOptions
	);
	CSLDestroy(papszOptions);

	if (!p_dataset) {
		throw std::runtime_error("unable to create dataset with driver.");
	}

	CPLErr err = p_dataset->SetGeoTransform(geotransform);
	if (err) {
		throw std::runtime_error("error setting geotransform.");
	}

	err = p_dataset->SetProjection(projection.c_str());
	if (err) {
		throw std::runtime_error("error setting projection.");
	}

	for (size_t band = 0; band < bandCount; band++) {
		GDALRasterBand *p_band = p_dataset->GetRasterBand(band + 1);
		p_band->SetDescription(bands[band].name.c_str());
		p_band->SetNoDataValue(bands[band].nan);
		p_band->GetBlockSize(&bands[band].xBlockSize, &bands[band].yBlockSize);
		bands[band].p_band = p_band;
	}

	return p_dataset;
}


/**
 *
 */
inline void
addBandToMEMDataset(
	GDALDataset *p_dataset,
	RasterBandMetaData& band)
{
	band.p_buffer = VSIMalloc3(
		p_dataset->GetRasterXSize(),
		p_dataset->GetRasterYSize(),
		band.size
	);

	CPLErr err;
	char **papszOptions = nullptr;
	const char *datapointer = std::to_string((size_t)band.p_buffer).c_str();
	papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", datapointer);
	
	err = p_dataset->AddBand(band.type, papszOptions);
	CSLDestroy(papszOptions);
	if (err) {
		throw std::runtime_error("unable to add band to dataset.");
	}
	
	GDALRasterBand *p_band = p_dataset->GetRasterBand(p_dataset->GetRasterCount());
	p_band->SetNoDataValue(band.nan);
	p_band->SetDescription(band.name.c_str());
	band.p_band = p_band;
}

/**
 *
 */
inline void
addBandToVRTDataset(
	GDALDataset *p_dataset,
	RasterBandMetaData& band,
	std::string tempFolder,
	std::string key,
	std::vector<VRTBandDatasetInfo>& VRTBandInfo)
{
	std::filesystem::path tmpPath = tempFolder;
	std::filesystem::path tmpName = "strat_breaks_" + key + ".tif";
	tmpPath = tmpPath / tmpName;

	VRTBandDatasetInfo info;
	info.filename = tmpPath.string();
	
	//there may be errors from GDAL if trying to tile a raster with scanline blocks
	//ensure these errors don't happen by not tiling when block sizes are a scanline
	bool useTiles = band.xBlockSize != p_dataset->GetRasterXSize() && 
					band.yBlockSize != p_dataset->GetRasterYSize();
	
	double geotransform[6];
	CPLErr err = p_dataset->GetGeoTransform(geotransform);
	if (err) {
		throw std::runtime_error("unable to get geotransform from dataset.");
	}

	//create the VRT band sub-dataset as a Geotiff
	info.p_dataset = createDataset(
		info.filename,
		"GTiff",
		p_dataset->GetRasterXSize(),
		p_dataset->GetRasterYSize(),
		geotransform,
		std::string(p_dataset->GetProjectionRef()),
		&band,
		1,
		useTiles
	);

	//add the sub-datasets band to the bands vector
	GDALRasterBand *p_band = info.p_dataset->GetRasterBand(1);
	p_band->SetDescription(band.name.c_str());
	p_band->SetNoDataValue(band.nan);
	p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
	band.p_band = p_band;
		
	VRTBandInfo.push_back(info);	
}


/**
 *
 */
inline void
addBandToVRTDataset(
	GDALDataset *p_dataset,
	RasterBandMetaData& band,
	VRTBandDatasetInfo& info
) {
	//adjust the options for the VRT band
	char **papszOptions = nullptr;

	//there may be errors from GDAL if trying to tile a raster with scanline blocks
	//ensure these errors don't happen by not tiling when block sizes are a scanline
	bool useTiles = band.xBlockSize != p_dataset->GetRasterXSize() && 
			band.yBlockSize != p_dataset->GetRasterYSize();

	if (useTiles) {
		const char *xBlockSize = std::to_string(band.xBlockSize).c_str();
		const char *yBlockSize = std::to_string(band.yBlockSize).c_str();
		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		papszOptions = CSLSetNameValue(papszOptions, "BLOCKXSIZE", xBlockSize);
		papszOptions = CSLSetNameValue(papszOptions, "BLOCKYSIZE", yBlockSize);
	}

	const char *filename = info.filename.c_str();
	papszOptions = CSLSetNameValue(papszOptions, "subclass", "VRTRawRasterBand");
	papszOptions = CSLSetNameValue(papszOptions, "SourceFilename", filename);

	//create the VRT band specifying the sub-dataset
	CPLErr err = p_dataset->AddBand(band.type, papszOptions);
	CSLDestroy(papszOptions);
	if (err) {
		throw std::runtime_error("unable to add band to dataset.");
	}

	GDALRasterBand *p_VRTBand = p_dataset->GetRasterBand(p_dataset->GetRasterCount());
	p_VRTBand->SetDescription(band.name.c_str());
	p_VRTBand->SetNoDataValue(band.nan);
}

/**
 *
 */
inline GDALDataset *
adjustBandsAndCreateDataset(
	std::string filename,
	std::string driver,
	int width,
	int height,
	double *geotransform,
	std::string projection,
	std::vector<RasterBandMetaData>& bands,
	size_t size,
	GDALDataType type,
	bool largeRaster)
{		
	//tiles must not be scanlines, as trying to set block size when they represent scanlines 
	//may result in GDAL errors due to the tile array being too large.
	bool useTiles = bands[0].xBlockSize != width &&
			bands[0].yBlockSize != height;

	for (size_t band = 0; band < bands.size(); band++) {
		bands[band].size = size;
		bands[band].type = type;
		bands[band].p_buffer = !largeRaster ? VSIMalloc3(height, width, size) : nullptr;
	}

	return createDataset(
		filename,
		driver,
		width,
		height,
		geotransform,
		projection,
		bands.data(),
		bands.size(),
		useTiles
	);
}

/**
 *
 */
inline void
rasterBandIO(
	RasterBandMetaData& band,
	void *p_buffer,
	int xBlockSize,
	int yBlockSize,
	int xBlock,
	int yBlock,
	int xValid,
	int yValid,
	bool read)
{
	CPLErr err;
	bool useBlock = xBlockSize == band.xBlockSize &&
						yBlockSize == band.yBlockSize;

	band.mutex.lock();
	if (useBlock && read) {
		err = read ?
			band.p_band->ReadBlock(xBlock, yBlock, p_buffer) :
			band.p_band->WriteBlock(xBlock, yBlock, p_buffer);
	}
	else {
		err = band.p_band->RasterIO(
			read ? GF_Read : GF_Write,
			xBlock * xBlockSize,
			yBlock * yBlockSize,
			xValid,
			yValid,
			p_buffer,
			xBlockSize,
			yBlockSize,
			band.type,
			0,
			0
		);
	}
	band.mutex.unlock();
	
	if (err) {
		throw read ?
			std::runtime_error("unable to read block from raster.") :
			std::runtime_error("unable to write block to raster.");
	}
}
