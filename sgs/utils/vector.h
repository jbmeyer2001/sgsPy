#pragma once

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <ogr_core.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

/**
 *
 */
class GDALVectorWrapper {
	private:
	GDALDatasetUniquePtr p_dataset;

	public:
	
	/**
	 *
	 */	
	GDALVectorWrapper(std::string filename);

	/**
	 *
	 */
	std::vector<std::string> getLayers();

	/**
	 *
	 */
	std::unordered_map<std::string, std::string> getLayerInfo(std::string layerName);

	/**
	 *
	 */
	OGRLayer *getLayer(std::string layerName);
};

PYBIND11_MODULE(vector, m) {
	py::class_<GDALVectorWrapper>(m, "GDALVectorWrapper")
		.def(py::init<std::string>())
		.def("get_layers", &GDALVectorWrapper::getLayers)
		.def("get_layer_info", &GDALVectorWrapper::getLayerInfo);
}
