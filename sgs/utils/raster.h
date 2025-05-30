#pragma once

#include <gdal_priv.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

/**
 *
 */
class GDALRasterWrapper {
	private:
	GDALDatasetUniquePtr p_dataset;
	double geotransform[6];

	public:
	/**
	 *
	 */
	GDALRasterWrapper(std::string filename);

	/**
	 *
	 */
	std::string getDriver();

	/**
	 *
	 */
	py::dict getCRS();

	/**
	 *
	 */
	int getWidth();

	/**
	 *
	 */
	int getHeight();

	/**
	 *
	 */
	int getLayers();

	/**
	 *
	 */
	double getXMax();

	/**
	 *
	 */
	double getXMin();

	/**
	 *
	 */
	double getYMax();

	/**
	 *
	 */
	double getYMin();

	/**
	 *
	 */
	double getPixelHeight();

	/**
	 *
	 */
	double getPixelWidth();

	/**
	 *
	 */
	std::vector<std::string> getBands();
		
	/*
	template <typename T>
	//TODO find some way to make this accessable by the python [] operator
	std::vector<T> getRasterAsNumpy() {
		//TODO add
		return std::vector<int>;
	}

	std::vector<T> getVirtualMemoryRaster() {
		//TODO add
		return std::vector<int>;
	}
	*/
};

PYBIND11_MODULE(raster, m) {
	py::class_<GDALRasterWrapper>(m, "GDALRasterWrapper")
		.def(py::init<std::string>())
		.def("get_driver", &GDALRasterWrapper::getDriver)
		.def("get_crs", &GDALRasterWrapper::getCRS)
		.def("get_height", &GDALRasterWrapper::getHeight)
		.def("get_width", &GDALRasterWrapper::getWidth)
		.def("get_layers", &GDALRasterWrapper::getLayers)
		.def("get_xmin", &GDALRasterWrapper::getXMin)
		.def("get_xmax", &GDALRasterWrapper::getXMax)
		.def("get_ymin", &GDALRasterWrapper::getYMin)
		.def("get_ymax", &GDALRasterWrapper::getYMax)
		.def("get_pixel_height", &GDALRasterWrapper::getPixelHeight)
		.def("get_pixel_width", &GDALRasterWrapper::getPixelWidth)
		.def("get_bands", &GDALRasterWrapper::getBands);
		//.def( something for getting the raster as a numpy array
}
