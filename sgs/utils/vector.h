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
    	GDALVectorWrapper(GDALDataset *p_dataset) {
		this->p_dataset = GDALDatasetUniquePtr(p_dataset);
	}

	/**
	 * Constructor for GDALVectorWrapper class. This constructor is meant to be used when importing
	 * data from another Python geospatial library, like geopandas. The geometries are stored in a python
	 * list as python bytes objects, and passed to this constructor through pybind11 as a vector of strings.
	 *
	 * Additionally, the projection, geometry type (of all geometries), and layer name are all used to create
	 * the GDALVectorWrapper object.
	 *
	 * @param std::vector<std::string> geometries
	 * @param std::string projection
	 * @param std::string geomtype
	 * @param std::string name
	 */
	GDALVectorWrapper(std::vector<std::string> geometries, std::string projection, std::string geomtype, std::string name) {
		GDALAllRegister();

		//set geometry type
		OGRwkbGeometryType type = wkbUnknown;
		if (geomtype == "point") {
			type = wkbPoint;
		}		
		else if (geomtype == "line") {
			type = wkbLineString;
		}
		else if (geomtype == "polygon") {
			type = wkbPolygon;
		}
		else {
			throw std::runtime_error("geomtype is not one of the accepted inputs: 'point', 'line', or 'polygon'.");
		}

		//set spatial reference
	       	this->srs.importFromWkt(projection.c_str());

		//create dataset and layer
		GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
		if (!p_driver) {
			throw std::runtime_error("unable to create dataset driver");
		}
		GDALDataset *p_dataset = p_driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
		if (!p_dataset) {
			throw std::runtime_error("unable to create dataset from driver.");
		}
		OGRLayer *p_layer = p_dataset->CreateLayer(name.c_str(), &srs, type, nullptr);
		if (!p_layer) {
			throw std::runtime_error("unable to create dataset layer.");
		}
		this->p_dataset = GDALDatasetUniquePtr(p_dataset);
	
		//add geometries to dataset
		for (size_t i = 0; i < geometries.size(); i++) {
			size_t size = geometries[i].length();
			void *p_bytes = (void *)geometries[i].c_str();
			OGRGeometry *p_geometry;

			OGRErr ogrerr = OGRGeometryFactory::createFromWkb(p_bytes, &srs, &p_geometry, size);
		       	if (!p_geometry) {
				throw std::runtime_error("unable to create geometry from wkb");
			}

			switch(type) {
				case OGRwkbGeometryType::wkbPoint:
					if (wkbFlatten(p_geometry->getGeometryType()) == wkbPoint) {
						OGRPoint *p_point = p_geometry->toPoint();
						addPoint(p_point, p_layer);	
					}
					else if (wkbFlatten(p_geometry->getGeometryType()) == wkbMultiPoint) {
						for (auto& p_point : *p_geometry->toMultiPoint()) {
							addPoint(p_point, p_layer);
						}
					}
					else {
						throw std::runtime_error("error mixing multiple geometry types in sgs GDALVectorWrapper object.");
					}
					break;
				case OGRwkbGeometryType::wkbLineString:
					if (wkbFlatten(p_geometry->getGeometryType()) == wkbLineString) {
						OGRLineString *p_lineString = p_geometry->toLineString();
						addLineString(p_lineString, p_layer);
					}
					else if (wkbFlatten(p_geometry->getGeometryType()) == wkbMultiLineString) {
						for (auto& p_lineString : *p_geometry->toMultiLineString()) {
							addLineString(p_lineString, p_layer);
						}
					}
					else {
						throw std::runtime_error("error mixing multiple geometry types in sgs GDALVectorWrapper object.");
					}
					break;
				case OGRwkbGeometryType::wkbPolygon:
					if (wkbFlatten(p_geometry->getGeometryType()) == wkbPolygon) {
						OGRPolygon *p_polygon = p_geometry->toPolygon();
						addPolygon(p_polygon, p_layer);
					}
					else if (wkbFlatten(p_geometry->getGeometryType()) == wkbMultiPolygon) {
						for (auto& p_polygon : *p_geometry->toMultiPolygon()) {
							addPolygon(p_polygon, p_layer);
						}
					}
					else {
						throw std::runtime_error("error mixing multiple geometry types in sgs GDALVectorWrapper object.");
					}
					break;
				default:
					throw std::runtime_error("should not be here!!! (bug).");
			}

			OGRGeometryFactory::destroyGeometry(p_geometry);
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
};
