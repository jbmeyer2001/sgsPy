/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for vector operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#pragma once

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <ogr_core.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

/**
 * Wrapper class for GDAL dataset containing a vector image.
 *
 * This class provides getter methods for important vector data and metadata,
 * as well as access to individual layers as OGRLayer * pointers.
 *
 * The Python side of the application only has access to metadata on each layer.
 *
 * The C++ side of the application has access to the getLayer() function.
 */
class GDALVectorWrapper {
	private:
	GDALDatasetUniquePtr p_dataset;

	public:	
	/**
	 * Constructor for GDALVectorWrapper class. this method registers
	 * drivers, and creates a GDALDataset object.
	 *
	 * @param filename as std::string
	 * @throws std::runtime_error if dataset is not initialized
	 */	
	GDALVectorWrapper(std::string filename);

	/**
	 * Getter method for vector layer names.
	 *
	 * @returns std::vector<std::string> an array of layer names as strings
	 */
	std::vector<std::string> getLayerNames();

	/**
	 * Getter method for layer parameters. Including:
	 * feature count
	 * field count
	 * geometry type
	 * xmin
	 * xmax
	 * ymin
	 * ymax
	 *
	 * @param std::string layer name
	 * @returns std::unordered_map<std::string, std::string> layer info converted to Python dict automatically
	 */
	std::unordered_map<std::string, std::string> getLayerInfo(std::string layerName);

	/**
	 * Getter method for a particular layer.
	 *
	 * @param std::string layer name
	 * @returns OGRLayer * pointer to a created instance of OGRLayer
	 */
	OGRLayer *getLayer(std::string layerName);

	/**
	 * Getter method for the layer Geometries. Note that every geometry 
	 * within the layer must be of type Point or MultiPoint.
	 *
	 * The points are stored and returned as a 2d array of doubles
	 * (std::vector<std::vector<double>>)
	 * The X coordinates are stored in index 0, Y in 1.
	 * 
	 * The array should be indexed like so:
	 * arr[0 if x or 1 if y][point_index]
	 *
	 * For example, to get the x and y coordinates of the 3rd point: 
	 * x = arr[0][2]; 
	 * y = arr[1][2]
	 *
	 * @param std::string layer name
	 * @returns std::vector<std::vector<double>> 2D array of doubles representing points
	 * @throws std::runtime_error if an incorrect geometry is encountered
	 */
	std::vector<std::vector<double>> getPoints(std::string layerName);

	/**
	 * Getter method for the layer Geometries. Note that every geometry
	 * within the layer must be of type LineString or MultiLineString.
	 *
	 * The LineStrings are stored and returned as a 3d array of doubles.
	 * (std::vector<std::vector<std::vector<double>>>)
	 * Indexing the outer array gives the Line, which is a 2d array of points.
	 * The X coordinates are stored in index 0, Y in 1.
	 *
	 * The array should be indexed like so:
	 * arr[line_index][0 if x or 1 if y][point_index]
	 *
	 * For example, to get the x and y coordinates of the 2nd point in the 5th line:
	 * x = arr[4][0][1]; 
	 * y = arr[4][1][1];
	 *
	 * @param std::string layer name
	 * @return std::vector<std::vector<std::vector<double>>> 3d array of doubles representing LineStrings.
	 * @throws std::runtime_error if an incorrect geometry is encountered
	 */
	std::vector<std::vector<std::vector<double>>> getLineStrings(std::string layerName);
};
