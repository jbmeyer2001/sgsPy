/******************************************************************************
 *
 * Project: sgs
 * Purpose: generate an access mask using vector and raster datasets
 * Author: Joseph Meyer
 * Date: October, 2025
 *
 ******************************************************************************/

#pragma once

#include "vector.h"
#include "helper.h"

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <ogr_core.h>

/**
 * This struct handles existing sample plot points. It has a constructor
 * which takes a GDALVectorWrapper, geotransform, and width of the raster
 * as parameters. 
 *
 * The points are stored within an unordered map. This
 * map can then be checked to see if it contains particular index values.
 */
struct Existing {
	bool used;
	std::unordered_map<int64_t, OGRPoint> samples;
	double IGT[6];
	int64_t width;

	/**
	 * Constructor for the Existing struct.
	 *
	 * First, check to ensure the GDALRasterWrapper isn't a null pointer.
	 *
	 * Calculates the inverse geotransform of the given transform, using 
	 * GDAL's inv_geotransform function. 
	 *
	 * Iterates through the features of the input vector. Every feature
	 * must be either a Point or a MultiPoint. Every point is converted
	 * from their x and y coordinates to a single index value, which
	 * is then stored in an unordered map. The existing points are
	 * also added to the output layer.
	 *
	 * The width of the raster is used to calculate this final index value.
	 *
	 * @param GDALVectorWrapper *p_vect
	 * @param double *GT
	 * @param int64_t width
	 * @param OGRLayer *p_samples
	 */
	Existing(
		GDALVectorWrapper *p_vect, 
		double *GT, 
		int64_t width, 
		OGRLayer *p_samples, 
		bool plot,
	       	std::vector<double>& xCoords,
		std::vector<double>& yCoords) 
	{
		if (!p_vect) {
			this->used = false;
			return;
		}
		
		this->width = width;

		std::vector<std::string> layerNames = p_vect->getLayerNames();
		if (layerNames.size() > 1) {
			throw std::runtime_error("the file containing existing sample points must have only a single layer.");
		}

		//invert geotransform so we can use IGT to convert from point to indexes
		GDALInvGeoTransform(GT, this->IGT);

		std::string name = layerNames[0];
		OGRLayer *p_layer = p_vect->getLayer(name);

		for (const auto& p_feature : *p_layer) {
			OGRGeometry *p_geometry = p_feature->GetGeometryRef();
			switch (wkbFlatten(p_geometry->getGeometryType())) {
				case OGRwkbGeometryType::wkbPoint: {
					OGRPoint *p_point = p_geometry->toPoint();
					int64_t index = point2index<int64_t>(p_point->getX(), p_point->getY(), IGT, width);
					this->samples.emplace(index, *p_point);
					if (p_samples) {
						addPoint(p_point, p_samples);
						if (plot) {
							xCoords.push_back(p_point->getX());
							yCoords.push_back(p_point->getY());
						}
					}
					break;
				}
				case OGRwkbGeometryType::wkbMultiPoint: {
					for (const auto& p_point : *p_geometry->toMultiPoint()) {
						int64_t index = point2index<int64_t>(p_point->getX(), p_point->getY(), IGT, width);
						this->samples.emplace(index, *p_point);
						if (p_samples) {
							addPoint(p_point, p_samples);
							if (plot) {
								xCoords.push_back(p_point->getX());
								yCoords.push_back(p_point->getY());
							}
						}
					}
					break;
				}
				default:
					throw std::runtime_error("the file containing existing sample points must have only Point or MultiPoint geometries.");
			}		
		}

		this->used = true;
	}

	/**
	 * Checker function which converts x and y indices values to a
	 * single index value. If this index is contained in the samples
	 * unordered_set, True is returned, otherwise the result will be 
	 * false.
	 *
	 * This function will be used when determining sample plot placement,
	 * when iterating through an input raster. 
	 */
	inline bool
	containsIndex(int64_t x, int64_t y) {
		int64_t index = y * this->width + x;
		return this->samples.find(index) != this->samples.end();	
	}

	/**
	 * if the map has already been checked, gets the OGRPoint which
	 * is associated with a particular value.
	 */
	inline OGRPoint
	getPoint(int64_t x, int64_t y) {
		int64_t index = y * this->width + x;
		return this->samples.find(index)->second;
	}

	/**
	 * Checker function which converts x coordinate and y coordinate
	 * values to an index usign the inverse geotransform.
	 *
	 * If this index is contained in the samples unordered_set, True
	 * is returned, otherwise the result will be false.
	 */
	inline bool
	containsCoordinates(double xCoord, double yCoord) {
		int64_t index = point2index<int64_t>(xCoord, yCoord, this->IGT, this->width);
		return this->samples.find(index) != this->samples.end();
	}

	/**
	 * Get the number of existing sample points.
	 */
	inline size_t
	count() {
		return this->samples.size();
	}
};
