/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALRasterWrapper pybind module
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

	
py::class_<sgs::raster::GDALRasterWrapper>(m, "GDALRasterWrapper")
	.def(py::init<std::string, std::string>())
	.def(py::init<py::buffer, std::vector<double>, std::string, std::vector<double>, std::vector<std::string>, std::string>())
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

