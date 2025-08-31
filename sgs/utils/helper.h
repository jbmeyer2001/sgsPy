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
#include <gdal_priv.h>

#define MAXINT8		127
#define MAXINT16	32767

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
createDataset(
	std::string filename,
	std::string driverName, 
	int width, 
	int height, 
	int bands,
	GDALDataType type,
	double *geotransform, 
	std::string projection) 
{
	GDALAllRegister();
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName(driverName.c_str());
	if (!p_driver) {
		throw std::runtime_error("unable to find dataset driver.");
	}

	GDALDataset *p_dataset = p_driver->Create(
		filename.c_str(),
		width,
		height,
		bands,
		type,
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
inline void
addBandToMEMDataset(
	std::vector<GDALRasterBand *>& bands,
	std::vector<void *>& stratBandBuffers,
	GDALDataset *p_dataset,
	GDALDataType type,
	size_t size,
	std::string bandName)
{
	void *p_strata = VSIMalloc3(
		p_dataset->GetRasterXSize(),
		p_dataset->GetRasterYSize(),
		size
	);
	stratBandBuffers.push_back(p_strata);

	CPLErr err;
	char **papszOptions = nullptr;
	papszOptions = CSLSetNameValue(
		papszOptions,
		"DATAPOINTER",
		std::to_string((size_t)p_strata).c_str()
	);

	err = p_dataset->AddBand(type, papszOptions);
	if (err) {
		throw std::runtime_error("unable to add band to dataset.");
	}

	bands.push_back(p_dataset->GetRasterBand(p_dataset->GetRasterCount()));
}

/**
 *
 */
inline void
addBandToVRTDataset(
	std::vector<GDALRasterBand *>& bands,
	std::vector<GDALDataset *>& VRTSubDatasets,
	GDALDataset * p_dataset,
	GDALDataType type,
	std::string bandName,
	int xBlockSize,
	int yBlockSize,
	std::string newFilename,
	double *geotransform,
	std::string projection)
{
	//set block size of new band
	CPLErr err;
	char **papszOptions = nullptr;
	papszOptions = CSLSetNameValue(
		papszOptions,
		"BLOCKXSIZE",
		std::to_string(xBlockSize).c_str()
	);
	papszOptions = CSLSetNameValue(
		papszOptions,
		"BLOCKYSIZE",
		std::to_string(yBlockSize).c_str()
	);

	//create the VRT band sub-dataset as a Geotiff
	GDALDataset *p_VRTSubDataset = createDataset(
		newFilename,
		"GTiff",
		p_dataset->GetRasterXSize(),
		p_dataset->GetRasterYSize(),
		1,
		type,
		geotransform,
		projection
	);
	VRTSubDatasets.push_back(p_VRTSubDataset);

	//add the sub-datasets band to the bands vector
	GDALRasterBand *p_band = p_VRTSubDataset->GetRasterBand(1);
	p_band->SetDescription(bandName.c_str());
	p_band->SetNoDataValue(-1);
	bands.push_back(p_band);

	//adjust the options for the VRT band
	papszOptions = CSLSetNameValue(
		papszOptions,
		"subclass",
		"VRTRawRasterBand"
	);
	papszOptions = CSLSetNameValue(
		papszOptions,
		"SourceFilename",
		newFilename.c_str()
	);
		
	//create the VRT band specifying the sub-dataset
	err = p_dataset->AddBand(type, papszOptions);
	if (err) {
		throw std::runtime_error("unable to add band to dataset.");
	}

	GDALRasterBand *p_VRTBand = p_dataset->GetRasterBand(p_dataset->GetRasterCount());
	p_VRTBand->SetDescription(bandName.c_str());
	p_VRTBand->SetNoDataValue(-1);
}
