/******************************************************************************
 *
 * Project: sgs
 * Purpose: Create Python bindings from C++ code
 * Author: Joseph Meyer
 * Date: November, 2025
 *
 ******************************************************************************/

#include "utils/raster.h"
#include "utils/vector.h"
#include "calculate/pca/pca.h"
#include "sample/clhs/clhs.h"
#include "sample/srs/srs.h"
#include "sample/strat/strat.h"
#include "sample/systematic/systematic.h"
#include "stratify/breaks/breaks.h"
#include "stratify/map/map.h"
#include "stratify/poly/poly.h"
#include "stratify/quantiles/quantiles.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(_sgs, m) {
	// source code in sgs/utils/raster.h
	py::class_<sgs::raster::GDALRasterWrapper>(m, "GDALRasterWrapper")
		.def(py::init<std::string>())
		.def(py::init<py::buffer, std::vector<double>, std::string, std::vector<double>, std::vector<std::string>>())
		.def("get_driver", &sgs::raster::GDALRasterWrapper::getDriver)
		.def("get_crs", &sgs::raster::GDALRasterWrapper::getCRS)
		.def("get_projection", &sgs::raster::GDALRasterWrapper::getFullProjectionInfo)
		.def("get_height", &sgs::raster::GDALRasterWrapper::getHeight)
		.def("get_width", &sgs::raster::GDALRasterWrapper::getWidth)
		.def("get_band_count", &sgs::raster::GDALRasterWrapper::getBandCount)
		.def("get_xmin", &sgs::raster::GDALRasterWrapper::getXMin)
		.def("get_xmax", &sgs::raster::GDALRasterWrapper::getXMax)
		.def("get_ymin", &sgs::raster::GDALRasterWrapper::getYMin)
		.def("get_ymax", &sgs::raster::GDALRasterWrapper::getYMax)
		.def("get_pixel_height", &sgs::raster::GDALRasterWrapper::getPixelHeight)
		.def("get_pixel_width", &sgs::raster::GDALRasterWrapper::getPixelWidth)
		.def("get_bands", &sgs::raster::GDALRasterWrapper::getBands)
		.def("get_band_nodata_value", &sgs::raster::GDALRasterWrapper::getBandNoDataValue)
		.def("get_raster_as_memoryview", &sgs::raster::GDALRasterWrapper::getRasterBandAsMemView)
		.def("get_raster_band_type_size", &sgs::raster::GDALRasterWrapper::getRasterBandTypeSize)
		.def("get_geotransform", &sgs::raster::GDALRasterWrapper::getGeotransformArray)
		.def("get_data_type", &sgs::raster::GDALRasterWrapper::getDataType)
		.def("set_temp_dir", &sgs::raster::GDALRasterWrapper::setTempDir)
		.def("get_temp_dir", &sgs::raster::GDALRasterWrapper::getTempDir)
		.def("release_band_buffers", &sgs::raster::GDALRasterWrapper::releaseBandBuffers)
		.def("close", &sgs::raster::GDALRasterWrapper::close);

	// source code in sgs/utils/vector.h
	py::class_<sgs::vector::GDALVectorWrapper>(m, "GDALVectorWrapper")
		.def(py::init<std::string>())
		.def(py::init<py::bytes, std::string, std::string>())
		.def("get_layer_names", &sgs::vector::GDALVectorWrapper::getLayerNames)
		.def("get_layer_info", &sgs::vector::GDALVectorWrapper::getLayerInfo)
		.def("get_points", &sgs::vector::GDALVectorWrapper::getPoints)
		.def("get_wkt_points", &sgs::vector::GDALVectorWrapper::getPointsAsWkt)
		.def("get_linestrings", &sgs::vector::GDALVectorWrapper::getLineStrings)
		.def("write", &sgs::vector::GDALVectorWrapper::write)
		.def("get_projection", &sgs::vector::GDALVectorWrapper::getFullProjectionInfo);

	// source code in sgs/calculate/pca/pca.h
	m.def("pca_cpp", &sgs::pca::pca);

	// source code in sgs/sample/clhs/clhs.h
	m.def("clhs_cpp", &sgs::clhs::clhs,
		pybind11::arg("p_raster"),
		pybind11::arg("nSamp"),
		pybind11::arg("iterations"),
		pybind11::arg("p_access"),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("tempFolder"),
		pybind11::arg("filename"));

	// source code in sgs/sample/srs/srs.h
	m.def("srs_cpp", &sgs::srs::srs, 
		pybind11::arg("p_raster"),
		pybind11::arg("numSamples"),
		pybind11::arg("mindist"),
		pybind11::arg("p_existing").none(true),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("tempFolder"),
		pybind11::arg("filename"));

	// source code in sgs/sample/strat/strat.h
	m.def("strat_cpp", &sgs::strat::strat,
		pybind11::arg("p_raster"),
		pybind11::arg("bandNum"),
		pybind11::arg("numSamples"),
		pybind11::arg("numStrata"),
		pybind11::arg("allocation"),
		pybind11::arg("weights"),
		pybind11::arg("p_mraster").none(true),
		pybind11::arg("mrastBandNum"),
		pybind11::arg("method"),
		pybind11::arg("wrow"),
		pybind11::arg("wcol"),
		pybind11::arg("mindist"),
		pybind11::arg("p_existing").none(true),
		pybind11::arg("force"),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("filename"),
		pybind11::arg("tempFolder"));

	// source code in sgs/sample/systematic/systematic.h
	m.def("systematic_cpp", &sgs::systematic::systematic,
		pybind11::arg("p_raster"),
		pybind11::arg("cellSize"),
		pybind11::arg("shape"),
		pybind11::arg("location"),
		pybind11::arg("p_existing").none(true),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("force"),
		pybind11::arg("plot"),
		pybind11::arg("filename"));

	// source code in sgs/stratify/breaks/breaks.h
	m.def("breaks_cpp", &sgs::breaks::breaks);

	// source code in sgs/stratify/map/map_stratifications.h
	m.def("map_cpp", &sgs::map::map);

	// source code in sgs/stratify/poly/poly.h
	m.def("poly_cpp", &sgs::poly::poly);

	// source code in sgs/stratify/quantiles/quantiles.h
	m.def("quantiles_cpp", &sgs::quantiles::quantiles);
}
