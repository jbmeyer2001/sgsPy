/******************************************************************************
 *
 * Project: sgs
 * Purpose: generate an access mask using vector and raster datasets
 * Author: Joseph Meyer
 * Date: July, 2025
 *
 ******************************************************************************/

#include "raster.h"
#include "vector.h"
#include "gdal_utils.h"

/*
 * Calculates a polygon for each LineString or MultiLineString in
 * the provided access polygon, Buffer() with outerBuffer, and with
 * innerBuffer if innerBuffer is not equal to zero.
 *
 * Takes the union of all outerBuffer polygons, and the union of all
 * innerBuffer polygons, then finds the difference. Get a single 
 * geometry (polygon) which is the difference.
 *
 * Rasterize the polygon into an access mask of type GDT_Byte,
 * return the dataset which contains this band.
 *
 * All geometries in the layer must be of type LineString or MultiLineString.
 *
 * @param std::string layerName layer to get access polygons for
 * @param double buffInner buffer not to sample from
 * @param double buffOuter buffer which must be sampled from
 * @returns GDALDataset * access mask
 */
GDALDataset *getAccessMask(
	GDALVectorWrapper *p_vector,
	GDALRasterWrapper *p_raster,
	std::string layerName, 
	double buffInner, 
	double buffOuter
); 
