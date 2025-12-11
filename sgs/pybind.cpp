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

PYBIND11_MODULE(_sgs, m) {
	// source code in sgs/utils/raster.h
	py::class_<GDALRasterWrapper>(m, "GDALRasterWrapper")
		.def(py::init<std::string>())
		.def(py::init<py::buffer, std::vector<double>, std::string, double, std::vector<std::string>>())
		.def("get_driver", &GDALRasterWrapper::getDriver)
		.def("get_crs", &GDALRasterWrapper::getCRS)
		.def("get_projection", &GDALRasterWrapper::getFullProjectionInfo)
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
		.def("get_raster_band_type_size", &GDALRasterWrapper::getRasterBandTypeSize)
		.def("get_geotransform", &GDALRasterWrapper::getGeotransformArray)
		.def("get_data_type", &GDALRasterWrapper::getDataType)
		.def("set_temp_dir", &GDALRasterWrapper::setTempDir)
		.def("get_temp_dir", &GDALRasterWrapper::getTempDir)
		.def("release_band_buffers", &GDALRasterWrapper::releaseBandBuffers)
		.def("close", &GDALRasterWrapper::close);

	// source code in sgs/utils/vector.h
	py::class_<GDALVectorWrapper>(m, "GDALVectorWrapper")
		.def(py::init<std::string>())
		.def("get_layer_names", &GDALVectorWrapper::getLayerNames)
		.def("get_layer_info", &GDALVectorWrapper::getLayerInfo)
		.def("get_points", &GDALVectorWrapper::getPoints)
		.def("get_wkt_points", &GDALVectorWrapper::getPointsAsWkt)
		.def("get_linestrings", &GDALVectorWrapper::getLineStrings);

	// source code in sgs/calculate/pca/pca.h
	m.def("pca_cpp", &pca::pca);

	// source code in sgs/sample/clhs/clhs.h
	m.def("clhs_cpp", &clhs::clhs,
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
	m.def("srs_cpp", &srs::srs, 
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
	m.def("strat_cpp", &strat::strat,
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
	m.def("systematic_cpp", &systematic::systematic,
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
	m.def("breaks_cpp", &breaks::breaks);

	// source code in sgs/stratify/map/map_stratifications.h
	m.def("map_cpp", &map::map);

	// source code in sgs/stratify/poly/poly.h
	m.def("poly_cpp", &poly::poly);

	// source code in sgs/stratify/quantiles/quantiles.h
	m.def("quantiles_cpp", &quantiles::quantiles);
}
