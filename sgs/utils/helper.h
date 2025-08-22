/******************************************************************************
 *
 * Project: sgs
 * Purpose: helper functions
 * Author: Joseph Meyer
 * Date: August, 2025
 *
 ******************************************************************************/

#pragma once

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
inline double
getPixelValueDependingOnType(
	GDALDataType type,
	void *p_data,
	size_t index
) {
	switch (type) {
		case GDT_Int8:
			return static_cast<double>(((int8_t *)p_data)[index]);
		case GDT_UInt16:
			return static_cast<double>(((uint16_t *)p_data)[index]);
		case GDT_Int16:
			return static_cast<double>(((int16_t *)p_data)[index]);
		case GDT_UInt32:
			return static_cast<double>(((uint32_t *)p_data)[index]);
		case GDT_Int32:
			return static_cast<double>(((int32_t *)p_data)[index]);
		case GDT_Float32:
			return static_cast<double>(((float *)p_data)[index]);
		case GDT_Float64:
			return ((double *)p_data)[index];
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
