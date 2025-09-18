/******************************************************************************
 *
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
			      GDALVectorWrapper()
******************************************************************************/
GDALVectorWrapper::GDALVectorWrapper(GDALDataset *p_dataset) {
	this->p_dataset = GDALDatasetUniquePtr(p_dataset);
}

/******************************************************************************
				  getDataset()
******************************************************************************/
GDALDataset *GDALVectorWrapper::getDataset() {
	return this->p_dataset.get();
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
std::unordered_map<std::string, std::string> 
GDALVectorWrapper::getLayerInfo(std::string layerName) {
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
	if (!p_layer->GetSpatialRef()) {
		std::cout << "WARNING: cannot get spatial reference for layer " << layerName << "." << std::endl;
	}
	else {
		retval.emplace("crs", std::string(p_layer->GetSpatialRef()->GetName()));
	}
	return retval;
}

/******************************************************************************
				   getLayer()
******************************************************************************/
OGRLayer *GDALVectorWrapper::getLayer(std::string layerName) {
	return this->p_dataset->GetLayerByName(layerName.c_str());
}



/******************************************************************************
				  getPoints()
******************************************************************************/
std::vector<std::vector<double>>
GDALVectorWrapper::getPoints(std::string layerName) {
	OGRLayer *p_layer = this->p_dataset->GetLayerByName(layerName.c_str());
	std::vector<double> xCoords;
	std::vector<double> yCoords;

	for (const auto& p_feature : *p_layer) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		switch (wkbFlatten(p_geometry->getGeometryType())) {
			case OGRwkbGeometryType::wkbPoint: {
				OGRPoint *p_point = p_geometry->toPoint();
				xCoords.push_back(p_point->getX());
				yCoords.push_back(p_point->getY());
				break;
			}
			case OGRwkbGeometryType::wkbMultiPoint: {
				for (const auto& p_point : *p_geometry->toMultiPoint()) {
					xCoords.push_back(p_point->getX());
					yCoords.push_back(p_point->getY());
				}
				break;
			}
			default:
				throw std::runtime_error("encountered a geometry which was not of type Point or MultiPoint.");
		}
	}

	return {xCoords, yCoords};
}

/******************************************************************************
				  getPointsAsWkt()
******************************************************************************/
std::vector<std::string>
GDALVectorWrapper::getPointsAsWkt(std::string layerName) {
	OGRLayer *p_layer = this->p_dataset->GetLayerByName(layerName.c_str());
	std::vector<std::string> retval;

	for (const auto& p_feature : *p_layer) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		switch (wkbFlatten(p_geometry->getGeometryType())) {
			case OGRwkbGeometryType::wkbPoint: {
				OGRPoint *p_point = p_geometry->toPoint();
				retval.push_back(p_point->exportToWkt());
				break;
			}
			case OGRwkbGeometryType::wkbMultiPoint: {
				for (const auto& p_point : *p_geometry->toMultiPoint()) {
					retval.push_back(p_point->exportToWkt());
				}
				break;
			}
			default:
				throw std::runtime_error("encountered a geometry which was not of type Point or MultiPoint.");
		}
	}

	return retval;
}

/******************************************************************************
				getLineStrings()
******************************************************************************/
std::vector<std::vector<std::vector<double>>> 
GDALVectorWrapper::getLineStrings(std::string layerName) {
	OGRLayer *p_layer = this->p_dataset->GetLayerByName(layerName.c_str());
	std::vector<std::vector<std::vector<double>>> retval;

	for (const auto& p_feature : *p_layer) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		switch (wkbFlatten(p_geometry->getGeometryType())) {
			case OGRwkbGeometryType::wkbLineString: {
				std::vector<double> xCoords;
				std::vector<double> yCoords;
				for (const auto& p_point : *p_geometry->toLineString()) {
					xCoords.push_back(p_point.getX());
					yCoords.push_back(p_point.getY());
				}
				retval.push_back({xCoords, yCoords});
				break;
			}
			case OGRwkbGeometryType::wkbMultiLineString: {
				for (const auto& p_lineString : *p_geometry->toMultiLineString()) {
					std::vector<double> xCoords;
					std::vector<double> yCoords;
					for (const auto& p_point : *p_lineString) {
						xCoords.push_back(p_point.getX());
						yCoords.push_back(p_point.getY());
					}
					retval.push_back({xCoords, yCoords});
				}
				break;
			}
			default:
				throw std::runtime_error("encountered a point which was not of type LineString or MultiLineString");
		}
	}

	return retval;
}

/******************************************************************************
				    write()
******************************************************************************/
void GDALVectorWrapper::write(std::string filename) {
	std::filesystem::path filepath = filename;
	std::string extension = filepath.extension().string();
		
	GDALAllRegister();
	GDALDriver *p_driver; 
	if (extension == ".geojson") {
		p_driver = GetGDALDriverManager()->GetDriverByName("GeoJSON");
	}
	else if (extension == ".shp") {
		p_driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	}
	else {
		throw std::runtime_error("file extension must be one of : .geojson, .shp");
	}

	GDALDataset *datasetCopy = p_driver->CreateCopy(
		filename.c_str(),
		this->p_dataset.get(),
		FALSE,
		nullptr,
		nullptr,
		nullptr
	);

	if (!datasetCopy) {
		std::cout << "failed to create dataset with filename " << filename << "." << std::endl;
	}

	CPLErr err = GDALClose(datasetCopy);

	if (err != CE_None) {
		std::cout << "failed to close dataset of file " << filename << ". The file output may not be correct." << std::endl;
	}
}
