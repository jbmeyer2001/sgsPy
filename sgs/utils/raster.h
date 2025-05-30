#pragma once

#include <gdal_priv.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

/**
 *
 */
class SpatialRaster {
	private:
	GDALDatasetUniquePtr p_dataset;
	double geotransform[6];

	public:
	/**
	 *
	 */
	SpatialRaster(std::string filename);

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
	py::class_<SpatialRaster>(m, "raster")
		.def(py::init<std::string>())
		.def("driver", &SpatialRaster::getDriver)
		.def("crs", &SpatialRaster::getCRS)
		.def("height", &SpatialRaster::getHeight)
		.def("width", &SpatialRaster::getWidth)
		.def("layers", &SpatialRaster::getLayers)
		.def("xmin", &SpatialRaster::getXMin)
		.def("xmax", &SpatialRaster::getXMax)
		.def("ymin", &SpatialRaster::getYMin)
		.def("ymax", &SpatialRaster::getYMax)
		.def("pixel_height", &SpatialRaster::getPixelHeight)
		.def("pixel_width", &SpatialRaster::getPixelWidth)
		.def("bands", &SpatialRaster::getBands);
		//.def( something for getting the raster as a numpy array
}
