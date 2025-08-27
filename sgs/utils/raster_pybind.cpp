/******************************************************************************
 *
 * Project: sgs
 * Purpose: Bind GDALDataset raster operations to Python
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"

/** 
 * pybind11 module for exposing GDALRasterWrapper data to Python.
 * Define only constructor and relevant getter functions.
 */
PYBIND11_MODULE(raster, m) {
	py::class_<GDALRasterWrapper>(m, "GDALRasterWrapper")
		.def(py::init<std::string>())
		.def("get_driver", &GDALRasterWrapper::getDriver)
		.def("get_crs", &GDALRasterWrapper::getCRS)
		.def("get_height", &GDALRasterWrapper::getHeight)
		.def("get_width", &GDALRasterWrapper::getWidth)
		.def("get_band_count", &GDALRasterWrapper::getBandCount)
		.def("get_xmin", &GDALRasterWrapper::getXMin)
		.def("get_xmax", &GDALRasterWrapper::getXMax)
		.def("get_ymin", &GDALRasterWrapper::getYMin)
		.def("get_ymax", &GDALRasterWrapper::getYMax)
		.def("get_pixel_height", &GDALRasterWrapper::getPixelHeight)
		.def("get_pixel_width", &GDALRasterWrapper::getPixelWidth)
		.def("get_bands", &GDALRasterWrapper::getBands)
		.def("get_band_nodata_value", &GDALRasterWrapper::getBandNoDataValue)
		.def("get_raster_as_memoryview", &GDALRasterWrapper::getRasterBandAsMemView)
		.def("get_raster_band_type_size", &GDALRasterWrapper::getRasterBandTypeSize);
}
