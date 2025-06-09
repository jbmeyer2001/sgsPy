/******************************************************************************
 *
 * Project: sgs
 * Purpose: Bind GDALDataset vector operations to Python
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "vector.h"

/**
 * pybind11 module for exposting GDALVectorWrapper data to Python.
 * Define only constructor and relevant getter functions.
 */
PYBIND11_MODULE(vector, m) {
	py::class_<GDALVectorWrapper>(m, "GDALVectorWrapper")
		.def(py::init<std::string>())
		.def("get_layer_names", &GDALVectorWrapper::getLayerNames)
		.def("get_layer_info", &GDALVectorWrapper::getLayerInfo)
		.def("get_geometries", &GDALVectorWrapper::getGeometries);
}
