/******************************************************************************
 *
 * Project: sgs
 * Purpose: generate an access mask using vector and raster datasets
 * Author: Joseph Meyer
 * Date: July, 2025
 *
 ******************************************************************************/

#include "access.h"

/******************************************************************************
				allocateRaster()
******************************************************************************/
GDALDataset *getAccessMask(
	GDALVectorWrapper *p_vector,
	GDALRasterWrapper *p_raster,
	std::string layerName, 
	double buffInner, 
	double buffOuter) 
{
	//step X:
	OGRMultiPolygon *buffInnerPolygons = new OGRMultiPolygon;
	OGRMultiPolygon *buffOuterPolygons = new OGRMultiPolygon;
	OGRGeometry *p_polygonMask;

	//step X:
	for (const auto& p_feature : *p_vector->getLayer(layerName.c_str())) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		OGRwkbGeometryType type = wkbFlatten(p_geometry->getGeometryType());

		if (type != wkbLineString && type != wkbMultiLineString) {
			throw std::runtime_error("all geometries in layer must be LineString or MultiLineString.");
		}

		if (buffInner == 0) {
			buffOuterPolygons->addGeometry(p_geometry->Buffer(buffOuter));
		}
		else {
			buffInnerPolygons->addGeometry(p_geometry->Buffer(buffInner));
			buffOuterPolygons->addGeometry(p_geometry->Buffer(buffOuter));
		}
	}	

	//step X: 
	if (buffInner == 0) {
		p_polygonMask = buffOuterPolygons->UnionCascaded();
		free(buffOuterPolygons);
	}
	else {
		OGRGeometry *buffOuterUnion = buffOuterPolygons->UnionCascaded();
		OGRGeometry *buffInnerUnion = buffInnerPolygons->UnionCascaded();
		p_polygonMask = buffOuterUnion->Difference(buffInnerUnion);
		free(buffOuterPolygons);
		free(buffInnerPolygons);
		free(buffOuterUnion);
		free(buffInnerUnion);
	}

	//TODO: do I really need to call GDALAllRegister() every time?
	GDALAllRegister();

	//TODO: error check output of these
	
	//step X: create new GDAL dataset to rasterize as access mask
	GDALDataset *p_accessPolygonDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create(
		"",
		0,
		0,
		0,
		GDT_Unknown,
		nullptr
	);

	//step X: create a layer of type Polygon in the new dataset
	OGRLayer *p_layer = p_accessPolygonDataset->CreateLayer("access", nullptr, wkbPolygon, nullptr);
	OGRFieldDefn field("index", OFTInteger);
	p_layer->CreateField(&field);
	
	//step X: add the access polygon to the new layer
	OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
	p_feature->SetField("index", 0);
	p_feature->SetGeometry(p_polygonMask);
	p_layer->CreateFeature(p_feature); //error handling here???
	OGRFeature::DestroyFeature(p_feature);

	//step X: get required info from raster
	double xRes = p_raster->getPixelWidth();
	double yRes = p_raster->getPixelHeight();
	double xMin = p_raster->getXMin();
	double xMax = p_raster->getXMax();
	double yMin = p_raster->getYMin();
	double yMax = p_raster->getYMax();

	//step X: generate options list for rasterization	
	char **argv;

	//specify the burn value for the polygon
	argv = CSLAddString(argv, "-burn");
	argv = CSLAddString(argv, std::to_string(1).c_str());

	//specify the layer
	argv = CSLAddString(argv, "-l");
	argv = CSLAddString(argv, "access");

	//specify the initialization values for the rest of the rastr
	argv = CSLAddString(argv, "-init");
	argv = CSLAddString(argv, std::to_string(0).c_str());

	//specify resolution
	argv = CSLAddString(argv, "-tr");
	argv = CSLAddString(argv, std::to_string(xRes).c_str());
	argv = CSLAddString(argv, std::to_string(yRes).c_str());

	//specify data type
	argv = CSLAddString(argv, "-ot");
	argv = CSLAddString(argv, "Byte");
	
	//specify extent
	argv = CSLAddString(argv, "-te");
	argv = CSLAddString(argv, std::to_string(xMin).c_str());
	argv = CSLAddString(argv, std::to_string(yMin).c_str());
	argv = CSLAddString(argv, std::to_string(xMax).c_str());
	argv = CSLAddString(argv, std::to_string(yMax).c_str());

	//specify the output format as in-memory
	argv = CSLAddString(argv, "-of");
	argv = CSLAddString(argv, "MEM");

	GDALRasterizeOptions *options = GDALRasterizeOptionsNew(argv, nullptr);

	//step X: rasterize vector creating in-memory dataset
	GDALDataset *p_accessRasterDataset = (GDALDataset *)GDALRasterize("",
		nullptr,
		p_accessPolygonDataset,
		options,
		nullptr
	);

	//step X: free dynamically allocated rasterization options
	GDALRasterizeOptionsFree(options);

	//step X: free no longer used polygon dataset
	free(p_accessPolygonDataset);

	//step X: return access mask
	return p_accessRasterDataset;
}
