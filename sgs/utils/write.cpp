/******************************************************************************
 *
 * Project: sgs
 * Purpose: Write sample points to vector file
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "write.h"

/******************************************************************************
				    write()
******************************************************************************/
void writeSamplePoints(std::vector<OGRPoint>& points, std::string filename) {
	std::filesystem::path filepath = filename;
	std::string extension = filepath.extension().string();

	//register format drivers
	GDALAllRegister();

	//create vector driver using file extension
	GDALDriver *p_driver;
	if (extension == ".geojson") {
		p_driver = GetGDALDriverManager()->GetDriverByName("GeoJSON");
	}
	else if (extension == ".shp") {
		p_driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	}
	else {
		throw std::runtime_error("file extension must be one of: .gdb, .geojson, .shp");
	}

	//ensure vector driver was created
	if (p_driver == nullptr) {
		throw std::runtime_error("driver for " + extension + " file type not available.");
	}

	//create and check dataset
	GDALDataset *p_dataset = p_driver->Create(filepath.string().c_str(), 0, 0, 0, GDT_Unknown, nullptr );
	if (p_dataset == nullptr) {
		throw std::runtime_error("could not create dataset with driver.");
	}

	//create and check layer
	OGRLayer *p_layer = p_dataset->CreateLayer("samples", nullptr, wkbPoint, nullptr);
	if (p_layer == nullptr) {
		throw std::runtime_error("could not create layer in datset.");
	}

	//create and check field
	OGRFieldDefn field("index", OFTInteger);
	if (extension == ".shp") {
		field.SetWidth(32);
	}
	if (p_layer->CreateField(&field) != OGRERR_NONE) {
		throw std::runtime_error("could not create name field.");
	}

	//add feature geometries to vector
	for (int i = 0; i < points.size(); i++) {
		OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
		p_feature->SetField("index", i);
		p_feature->SetGeometry(&points[i]);
		if (p_layer->CreateFeature(p_feature) != OGRERR_NONE) {
			throw std::runtime_error("Failed to create feature in file.");
		}

		//clean up memory
		OGRFeature::DestroyFeature(p_feature);
	}

	//close to ensure everything is written correctly
	GDALClose(p_dataset);
}
