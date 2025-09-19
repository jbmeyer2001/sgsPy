/******************************************************************************
 *
 * Project: sgs
 * Purpose: helper functions
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#pragma once

#include <iostream>
#include <filesystem>
#include <mutex>
#include <gdal_priv.h>

#define MAXINT8		127
#define MAXINT16	32767

/**
 * The RasterBandMetaData struct stores information
 * on a particular raster band. It stores:
 *
 * GDALRasterBand *p_band: 
 * 	a pointer to the associated GDALRasterBand.
 * void *p_buffer:
 * 	a pointer to the whole raster band as a buffer,
 * 	if the raster is a size where storing the whole band
 * 	is possible.
 * GDALDataType:
 * 	the pixel type of the raster band.
 * size_t size:
 * 	the size of the pixel type of the raster band.
 * std::string name:
 * 	the band name.
 * double nan:
 *	the bands particular no data value.
 * int xBlockSize:
 * 	the x component of the raster bands block size.
 * int yBlockSize:
 * 	the y component of hte raster bands block size.
 * std::mutex *p_mutex:
 * 	a pointer to the mutex for the GDAL dataset
 * 	corresponding to the raster band. This must
 * 	be locked/unlocked when reading/writing the band.
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
	std::mutex *p_mutex = nullptr;
};

/**
 * When adding a raster band to a VRT dataset, the dataset
 * (which will be added as a band) must be fully formed.
 * VRT datasets are used as a virtual type to retain processed
 * raster data on disk for use by the package, if the
 * whole band won't fit into memory. As a result, the
 * dataset which is added to the VRT is not fully
 * formed when it is created, because the data which
 * will be written to it has not yet been processed.
 *
 * This struct holds the information of a dataset (usually tif) which
 * will be added as a band to a VRT dataset once it
 * has been populated with processed data. The necessary
 * information to retain on that dataset are:
 *
 * GDALDataset *p_dataset: 
 * 	the pointer to the GDALDataset object representing the dataset.
 * std::string filename:
 * 	the file name of the dataset, typically within a temporary
 * 	folder on disk.
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
 * Helper function which determines the smallest 
 * signed integer type and it's corresponding size
 * which can fit the maxStrata without overflow
 * errors.
 *
 * @param size_t maxStrata
 * @param GDALDataType *p_dtye
 * @param size_t *p_size
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
 * This helper function is used by stratification functions to create
 * a virtual dataset to add bands to. The virtual dataset will either
 * be a MEM dataset or in-memory dataset, or a VRT dataset. The 
 * VRT dataset will be chosen if the raster is large enough
 * that putting the whole raster into memory should not be attempted.
 *
 * The function begins by getting a GDALDriver from the driver name
 * (which must be one of "MEM" or "VRT").
 *
 * Then, a new GDALDataset is created with the width and height
 * required. Type, band number, filename, and other options are
 * not used because the bands will be dynamically added with
 * potentially different types per band, and this dataset
 * will not be written directly to disk.
 *
 * The geotransform and projection are set.
 *
 * @param std::string driverName
 * @param int width
 * @param int height
 * @param double *geotransform
 * @Param std::string projection
 */
inline GDALDataset *
createVirtualDataset(
	std::string driverName, 
	int width, 
	int height, 
	double *geotransform, 
	std::string projection) 
{
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
 * This helper function is used when the user provides a specific
 * filename, and creates a dataset using the driver associated with
 * the filename extension.
 *
 * Note, the driverName must be determined from the filename before
 * this function is called.
 *
 * First, a GDALDriver is created with the driver according to the
 * driverName parameter. Then, if tiles should be used, those
 * CSL options are added to a CSL options list which will be passed
 * when creating the Dataset. User-given driver options are also
 * added to teh CSL options list.
 *
 * A dataset is created using the filename, height, width, band count,
 * and options specified. Note, the type of every raster in the dataset
 * must be the same (this is true right now, as the only driver allowed
 * is tiff, this may be changed in the future as more drivers are added.
 *
 * The geotransform and projection are set.
 *
 * Finally, since the band information (specifically type) must be known
 * before a non-virtual dataset is created, the RasterBandMetaData
 * structs associated with each band have not been fully populated. 
 * Bands could not habve been allocated, created, etc. until after 
 * all of them had been iterated through. Now that we are able to 
 * create the dataset with its associated bands, we must now populate
 * the RasterBandMetaData structs with the remaining meta data
 * including block sizes and GDALRasterBand pointers. The name
 * and nodata value of each band is also updated.
 *
 * @param std::string filename
 * @param std::string driverName
 * @param int width
 * @param int height
 * @param double *geotransform
 * @param std::string projection,k
 * @param RasterBandMetaData *bands
 * @param size_t bandCount
 * @param bool useTiles
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
	bool useTiles,
	std::map<std::string, std::string>& driverOptions) 
{
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName(driverName.c_str());
	if (!p_driver) {
		throw std::runtime_error("unable to find dataset driver.");
	}

	char **papszOptions = nullptr;
	if (useTiles) {
		if (driverOptions.find("TILED") == driverOptions.end()) {
			papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		}
		if (driverOptions.find("BLOCKXSIZE") == driverOptions.end()) {
			const char *xBlockSizeOption = std::to_string(bands[0].xBlockSize).c_str();
			papszOptions = CSLSetNameValue(papszOptions, "BLOCKXSIZE", xBlockSizeOption);
		}
		if (driverOptions.find("BLOCKYSIZE") == driverOptions.end()) {
			const char *yBlockSizeOption = std::to_string(bands[0].yBlockSize).c_str();
			papszOptions = CSLSetNameValue(papszOptions, "BLOCKYSIZE", yBlockSizeOption);
		}
	}

	for (auto const& [key, val] : driverOptions) {
		papszOptions = CSLSetNameValue(papszOptions, key.c_str(), val.c_str());
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
 * This helper function dynamically adds a raster band to an existing
 * in-memory dataset (MEM).
 *
 * First, the buffer is allocated according to the height,
 * width, and pixel size of the entire band. This buffer
 * is passed as the 'datapointer' option to the AddBand
 * function so that the memory is now associated with
 * the dataset.
 *
 * Once the band has been created, it's name and nodata values
 * are set, and the RasterBandMetaData struct of the particular
 * raster band which was added is updated with a pointer to the
 * GDALRasterBand object.
 *
 * @param GDALDataset *p_dataset
 * @param RasterBandMetaData& band
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
 * This helper function is used to create a dataset which will
 * be used as the band of a VRT dataset. In order to pass the
 * checks for the addition of a dataset as a band to a VRT dataset
 * the file must be populated (written to). This function is used
 * before any data processing is undertaken to create an empty
 * tif file, so that data has a place to be written to, and is
 * added as a band to a VRT dataset once the processing has completed.
 *
 * When this function has completed there will be a new entry in
 * the VRTBandInfo vector containing a dataset pointer and filename.
 * The RasterBandMetaData object will also be updated with the band
 * information for the tif file created. It's worth noting that
 * the tmpPath is a temporary folder created for the SpatialRaster
 * which this function was called, that holds the .tif files which
 * make up the VRT dataset.
 *
 * first, the full file path is determined using the temp path,
 * as well as a unique key differentiating the different bands
 * which may be in the dataset. The useTiles boolean determines
 * whether the TILED, XBLOCKSIZE, and YBLOCKSIZE options
 * will be used in dataset creation. The reason why they wouldn't
 * be used is if the block size is a scanline -- this is because
 * GDAL will throw an error if the tile array would be larger
 * than 2GB, which would happen if tiles were scanlines on
 * some large images. Then, the geotransform is determined, 
 * and createDataset() is called.
 *
 * Finally, the new tif datasets only band is updated with
 * a name, and nodata value, and the corresponding 
 * RasterBandMetaData object is updated.
 *
 * @param GDALDataset *p_dataset
 * @param RasterBandMetaData& data
 * @param std::string tempFolder
 * @param std::string key
 * @param std::vector<VRTBandDatasetInfo>& VRTBandInfo
 */
inline void
createVRTBandDataset(
	GDALDataset *p_dataset,
	RasterBandMetaData& band,
	std::string tempFolder,
	std::string key,
	std::vector<VRTBandDatasetInfo>& VRTBandInfo,
	std::map<std::string, std::string>& driverOptions)
{
	std::filesystem::path tmpPath = tempFolder;
	std::filesystem::path tmpName = "strat_breaks_" + key + ".tif";
	tmpPath = tmpPath / tmpName;

	VRTBandDatasetInfo info;
	info.filename = tmpPath.string();
	
	bool useTiles = band.xBlockSize != p_dataset->GetRasterXSize() && 
			band.yBlockSize != p_dataset->GetRasterYSize();

	double geotransform[6];
	CPLErr err = p_dataset->GetGeoTransform(geotransform);
	if (err) {
		throw std::runtime_error("unable to get geotransform from dataset.");
	}

	info.p_dataset = createDataset(
		info.filename,
		"GTiff",
		p_dataset->GetRasterXSize(),
		p_dataset->GetRasterYSize(),
		geotransform,
		std::string(p_dataset->GetProjectionRef()),
		&band,
		1,
		useTiles,
		driverOptions
	);

	GDALRasterBand *p_band = info.p_dataset->GetRasterBand(1);
	p_band->SetDescription(band.name.c_str());
	p_band->SetNoDataValue(band.nan);
	p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
	band.p_band = p_band;
		
	VRTBandInfo.push_back(info);	
}


/**
 * This helper function adds an existing raster dataset as a band
 * to an existing VRT dataset.
 *
 * First, creation options are determined. The filename, as well as the VRT band
 * subclass are added as options. Then, the band is added and updated with 
 * a name and nodata value.
 *
 * @param GDALDataset *p_dataset
 * @param RasterBandMetaData& band
 * @param VRTBandDatasetInfo& info
 */
inline void
addBandToVRTDataset(
	GDALDataset *p_dataset,
	RasterBandMetaData& band,
	VRTBandDatasetInfo& info
) {
	char **papszOptions = nullptr;
	const char *filename = info.filename.c_str();
	papszOptions = CSLSetNameValue(papszOptions, "subclass", "VRTRawRasterBand");
	papszOptions = CSLSetNameValue(papszOptions, "SourceFilename", filename);

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
 * This helper function is used to either read or write a chunk
 * of memory from or to a raster band. If the block size of
 * the band corresponds to the block size of the memory,
 * ReadBlock() or WriteBlock() may be used. Otherwise,
 * RasterIO is used. The boolean parameter 'read' lets
 * the function know whether the band should be read
 * to the buffer, or the buffer should be written to the
 * band.
 *
 * @param RasterBandMetaData& band
 * @param void *p_buffer
 * @param int xBlockSize
 * @param int yBlockSize
 * @param int xBlock
 * @param int yBlock
 * @param int xValid
 * @param int yValid
 * @param bool read
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

	band.p_mutex->lock();
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
	band.p_mutex->unlock();
	
	if (err) {
		throw read ?
			std::runtime_error("unable to read block from raster.") :
			std::runtime_error("unable to write block to raster.");
	}
}
