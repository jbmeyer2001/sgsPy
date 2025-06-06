/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for vector operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "vector.h"

/******************************************************************************
			      GDALVectorWrapper()
******************************************************************************/
GDALVectorWrapper::GDALVectorWrapper(std::string filename) {
	//must register drivers before trying to open a dataset
	GDALAllRegister();

	//dataset
	this->p_dataset = GDALDatasetUniquePtr(GDALDataset::Open(filename.c_str(), GDAL_OF_VECTOR));
	if (!this->p_dataset) {
		throw std::runtime_error("dataset pointer is null after initialization, dataset unable to be initialized.");
	}
}

/******************************************************************************
				getLayerNames()
******************************************************************************/
std::vector<std::string> GDALVectorWrapper::getLayerNames() {
	std::vector<std::string> retval;

	for (OGRLayer *p_layer : this->p_dataset->GetLayers()) {
		retval.push_back(std::string(p_layer->GetName()));
	}

	return retval;
}

/******************************************************************************
				 getLayerInfo()
******************************************************************************/
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

/******************************************************************************
				   getLayer()
******************************************************************************/
OGRLayer *GDALVectorWrapper::getLayer(std::string layerName) {
	return this->p_dataset->GetLayerByName(layerName.c_str());
}
