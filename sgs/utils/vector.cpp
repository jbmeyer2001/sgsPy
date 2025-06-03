#include "vector.h"
#include <iostream>

GDALVectorWrapper::GDALVectorWrapper(std::string filename) {
	//must register drivers first
	GDALAllRegister();

	//dataset
	this->p_dataset = GDALDatasetUniquePtr(GDALDataset::Open(filename.c_str(), GDAL_OF_VECTOR));
}

std::vector<std::string> GDALVectorWrapper::getLayers() {
	std::vector<std::string> retval;

	for (OGRLayer *p_layer : this->p_dataset->GetLayers()) {
		retval.push_back(std::string(p_layer->GetName()));
	}

	return retval;
}

std::unordered_map<std::string, std::string> GDALVectorWrapper::getLayerInfo(std::string layerName) {
	std::unordered_map<std::string, std::string> retval;

	OGRLayer *p_layer = this->p_dataset->GetLayerByName(layerName.c_str());
	std::unique_ptr<OGREnvelope> extent = std::unique_ptr<OGREnvelope>(new OGREnvelope);
	p_layer->GetExtent(extent.get());

	retval.emplace("feature_count", std::to_string(p_layer->GetFeatureCount()));
	retval.emplace("field_count", std::to_string(p_layer->GetLayerDefn()->GetFieldCount()));
	retval.emplace("geometry_type", OGRGeometryTypeToName(p_layer->GetGeomType()));
	retval.emplace("xmin", std::to_string(extent->MinX));
	retval.emplace("xmax", std::to_string(extent->MaxX));
	retval.emplace("ymin", std::to_string(extent->MinY));
	retval.emplace("ymax", std::to_string(extent->MaxY));
	
	return retval;
}

OGRLayer *GDALVectorWrapper::getLayer(std::string layerName) {
	return this->p_dataset->GetLayerByName(layerName.c_str());
}
