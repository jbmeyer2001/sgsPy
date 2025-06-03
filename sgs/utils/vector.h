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
};

/**
 * pybind11 module for exposting GDALVectorWrapper data to Python.
 * Define only constructor and relevant getter functions.
 */
PYBIND11_MODULE(vector, m) {
	py::class_<GDALVectorWrapper>(m, "GDALVectorWrapper")
		.def(py::init<std::string>())
		.def("get_layers", &GDALVectorWrapper::getLayerNames)
		.def("get_layer_info", &GDALVectorWrapper::getLayerInfo);
}
