/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALVectorWrapper pybind module
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/
	
py::class_<sgs::vector::GDALVectorWrapper>(m, "GDALVectorWrapper")
	.def(py::init<std::string, std::string>())
	.def(py::init<py::bytes, std::string, std::string, std::string>())
	.def("get_layer_names", &sgs::vector::GDALVectorWrapper::getLayerNames)
	.def("get_layer_info", &sgs::vector::GDALVectorWrapper::getLayerInfo)
	.def("get_points", &sgs::vector::GDALVectorWrapper::getPoints)
	.def("get_wkt_points", &sgs::vector::GDALVectorWrapper::getPointsAsWkt)
	.def("get_linestrings", &sgs::vector::GDALVectorWrapper::getLineStrings)
	.def("write", &sgs::vector::GDALVectorWrapper::write)
	.def("get_projection", &sgs::vector::GDALVectorWrapper::getFullProjectionInfo);
