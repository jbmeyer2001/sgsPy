/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for vector operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#pragma once

#include <iostream>
#include <filesystem>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <ogr_core.h>
#include <ogr_geometry.h>
#include <ogr_api.h>
#include <ogr_recordbatch.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

/**
 * Wrapper class for GDAL dataset containing a vector image.
 *
 * This class provides getter methods for important vector data and metadata,
 * as well as access to individual layers as OGRLayer * pointers.
 *
 * The Python side of the application only has access to metadata on each layer.
 *
 * The C++ side of the application has access to the getLayer() function.
 */
class GDALVectorWrapper {
	private:
	GDALDatasetUniquePtr p_dataset;
	OGRSpatialReference srs; 
	bool haveSRS = false;

	public:	
	/**
	 * Constructor for GDALVectorWrapper class. this method registers
	 * drivers, and creates a GDALDataset object.
	 *
	 * @param filename as std::string
	 * @throws std::runtime_error if dataset is not initialized
	 */	
	GDALVectorWrapper(std::string filename) {
		//must register drivers before trying to open a dataset
		GDALAllRegister();

		//dataset
		this->p_dataset = GDALDatasetUniquePtr(GDALDataset::Open(filename.c_str(), GDAL_OF_VECTOR));
		if (!this->p_dataset) {
			throw std::runtime_error("dataset pointer is null after initialization, dataset unable to be initialized.");
		}
	}

	/**
	* Constructor for GDALVectorWrapper class.
	*
	* @param GDALDataset * pointer to existing GDAL dataset
	*/	
    	GDALVectorWrapper(GDALDataset *p_dataset, std::string projection) {
		this->p_dataset = GDALDatasetUniquePtr(p_dataset);
		OGRErr err = this->srs.importFromWkt(projection.c_str());
		if (err) {
			throw std::runtime_error("unable to get Spatial Reference System from projection string.");
		}
		this->haveSRS = true;
	}

	/**
	 * Constructor for GDALVectorWrapper class. This constructor is meant to be used when importing
	 * data from another Python geospatial library, like geopandas. The geodataset is converted
	 * to a geojson string and passed to the GDALOpenEx() function to create a GDALDataset. This
	 * dataset unfortunately won't have the correcct layer name, and may have an incorrect spatial
	 * reference system. Because of this, a GDALDataset is created with a new OGRLayer containing
	 * the correct layer name and spatial reference system. The geometries with their fields are
	 * then copied from the geojson-created dataset to the new initialized dataset. 
	 *
	 * @param std::vector<std::string> geometries
	 * @param std::string projection
	 * @param std::string name
	 */
	GDALVectorWrapper(std::string bytes, std::string projection, std::string name) {
		GDALAllRegister();
		
		//create dataset from geojson string to copy to new layer
		GDALDataset *p_indataset = GDALDataset::FromHandle(GDALOpenEx(bytes.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
		OGRLayer *p_inlayer = p_indataset->GetLayerByName("OGRGeoJSON");
	
		//set spatial reference
	       	OGRErr err = this->srs.importFromWkt(projection.c_str());
		if (err) {
			throw std::runtime_error("unable to get Spatial Reference System from projection string.");
		}
		this->haveSRS = true;

		//create dataset and layer with correct spatial reference and 
		GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
		if (!p_driver) {
			throw std::runtime_error("unable to create dataset driver.");
		}
		GDALDataset *p_dataset = p_driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
		if (!p_dataset) {
			throw std::runtime_error("unable to create dataset from driver.");
		}
		OGRLayer *p_outlayer = p_dataset->CreateLayer(name.c_str(), &this->srs, wkbUnknown, nullptr);
		if (!p_outlayer) {
			throw std::runtime_error("unable to create dataset layer.");
		}
		this->p_dataset = GDALDatasetUniquePtr(p_dataset);


		/* copy data over from dataset created using geojson string
		 *
		 * This is done, because when converting a geopandas geodataframe to json, the layer name,
		 * and (more importantly) the spatial reference system may not be transferred correctly.
		 * This way does include A LOT of copying data, however it ensures the srs is correct and
		 * all of the data from the geodataframe is moved over effectively.
		 */

		//copy over field definitions	
		OGRFeatureDefn *p_featdef = p_inlayer->GetLayerDefn();
		int fcount = p_featdef->GetFieldCount();
		for (int i = 0; i < fcount; i++) {
			OGRFieldDefn *p_fielddef = p_featdef->GetFieldDefn(i);
			OGRErr err = p_outlayer->CreateField(p_fielddef);
			if (err) {
				std::cout << "unable to copy field definition." << std::endl;
			}
		}

		//copy over features
		for (auto& p_infeature : p_inlayer) {
			OGRFeature *p_outfeature = OGRFeature::CreateFeature(p_featdef);
			for (int i = 0; i < fcount; i++) {
				p_outfeature->SetField(i, p_infeature->GetRawFieldRef(i));
			}
			OGRGeometry *p_geom = p_infeature->GetGeometryRef();
			p_geom->assignSpatialReference(&this->srs);
			p_outfeature->SetGeometry(p_geom);
			p_outlayer->CreateFeature(p_outfeature);
			OGRFeature::DestroyFeature(p_outfeature);
		}
	}	
			
	/**
	 * Getter method for the dataset pointer.
	 *
	 * @returns pointer to the underlying GDAL dataset
	 */
	GDALDataset *getDataset() {
		return this->p_dataset.get();
	}

	/**
	 * Getter method for vector layer names.
	 *
	 * @returns std::vector<std::string> an array of layer names as strings
	 */
	std::vector<std::string> getLayerNames() {
		std::vector<std::string> retval;

		for (OGRLayer *p_layer : this->p_dataset->GetLayers()) {
			retval.push_back(std::string(p_layer->GetName()));
		}

		return retval;
	}

	/**
	 * Getter method for layer parameters. Including:
	 * feature count
	 * field count
	 * geometry type
	 * xmin
	 * xmax
	 * ymin
	 * ymax
	 *
	 * @param std::string layer name
	 * @returns std::unordered_map<std::string, std::string> layer info converted to Python dict automatically
	 */
	std::unordered_map<std::string, std::string> getLayerInfo(std::string layerName) {
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

	/**
	 * Getter method for a particular layer.
	 *
	 * @param std::string layer name
	 * @returns OGRLayer * pointer to a created instance of OGRLayer
	 */
	OGRLayer *getLayer(std::string layerName) {
		return this->p_dataset->GetLayerByName(layerName.c_str());
	}

	/**
	 * Getter method for the layer Geometries. Note that every geometry 
	 * within the layer must be of type Point or MultiPoint.
	 *
	 * The points are stored and returned as a 2d array of doubles
	 * (std::vector<std::vector<double>>)
	 * The X coordinates are stored in index 0, Y in 1.
	 * 
	 * The array should be indexed like so:
	 * arr[0 if x or 1 if y][point_index]
	 *
	 * For example, to get the x and y coordinates of the 3rd point: 
	 * x = arr[0][2]; 
	 * y = arr[1][2]
	 *
	 * @param std::string layer name
	 * @returns std::vector<std::vector<double>> 2D array of doubles representing points
	 * @throws std::runtime_error if an incorrect geometry is encountered
	 */
	std::vector<std::vector<double>> getPoints(std::string layerName) {
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
	/**
	 * Getter method for the layer geometries, where every geometry
	 * must be of type Point or MultiPoint.
	 *
	 * the points are stored as wkt (well known text) strings in the
	 * return vector.
	 *
	 * @param std::string layer name
	 * @returns std::vector<std::string> vector of points as wkt
	 */
	std::vector<std::string> getPointsAsWkt(std::string layerName){
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

	/**
	 * Getter method for the layer Geometries. Note that every geometry
	 * within the layer must be of type LineString or MultiLineString.
	 *
	 * The LineStrings are stored and returned as a 3d array of doubles.
	 * (std::vector<std::vector<std::vector<double>>>)
	 * Indexing the outer array gives the Line, which is a 2d array of points.
	 * The X coordinates are stored in index 0, Y in 1.
	 *
	 * The array should be indexed like so:
	 * arr[line_index][0 if x or 1 if y][point_index]
	 *
	 * For example, to get the x and y coordinates of the 2nd point in the 5th line:
	 * x = arr[4][0][1]; 
	 * y = arr[4][1][1];
	 *
	 * @param std::string layer name
	 * @return std::vector<std::vector<std::vector<double>>> 3d array of doubles representing LineStrings.
	 * @throws std::runtime_error if an incorrect geometry is encountered
	 */
	std::vector<std::vector<std::vector<double>>> getLineStrings(std::string layerName){
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

	/**
	 * method for writing the existing dataset to a file.
	 *
	 * This is intended to be used to write in-memory datasets
	 * generated by sampling methods to disk, however it should 
	 * work for any dataset.
	 *
	 * @param std::string filename to write to
	 */
	void write(std::string filename) {
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

	/**
	 * Getter method for the full projection information as wkt.
	 *
	 * @returns std::string
	 */
	std::string getFullProjectionInfo() {
		char *p_proj;
		if (haveSRS) {
			srs.exportToPrettyWkt(&p_proj);
		}
		else {
			this->p_dataset->GetLayer(0)->GetSpatialRef()->exportToPrettyWkt(&p_proj);
		}
		
		std::string retval = std::string(p_proj);
		CPLFree(p_proj);
		return retval;
	}

	/**
	 * Getter method for OGRSpatialReference.
	 *
	 * @return OGRSpatialReference *
	 */
	OGRSpatialReference *getSRS(void) {
		if (!this->haveSRS) {
			throw std::runtime_error("do not have OGRSpatialReference associated with GDALVectorWrapper.");
		}

		return &this->srs;
	}
};
