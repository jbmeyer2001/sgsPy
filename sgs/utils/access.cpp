/******************************************************************************
 *
 * Project: sgs
 * Purpose: generate an access mask using vector and raster datasets
 * Author: Joseph Meyer
 * Date: July, 2025
 *
 ******************************************************************************/

#include "access.h"
#include <ogrsf_frmts.h>

/******************************************************************************
				getAccessMask()
******************************************************************************/
std::pair<GDALDataset *, void *>
getAccessMask(
	GDALVectorWrapper *p_vector,
	GDALRasterWrapper *p_raster,
	std::string layerName, 
	double buffInner, 
	double buffOuter) 
{
	//step 1: create multipolygon buffers
	OGRMultiPolygon *buffInnerPolygons = new OGRMultiPolygon;
	OGRMultiPolygon *buffOuterPolygons = new OGRMultiPolygon;
	OGRGeometry *p_polygonMask;

	//step 2: add geometries to access polygon buffers
	for (const auto& p_feature : *p_vector->getLayer(layerName.c_str())) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		OGRwkbGeometryType type = wkbFlatten(p_geometry->getGeometryType());

		switch (type) {
			case OGRwkbGeometryType::wkbLineString: {
				buffOuterPolygons->addGeometry(p_geometry->Buffer(buffOuter));
				if (buffInner != 0) {
					buffInnerPolygons->addGeometry(p_geometry->Buffer(buffInner));
				}
				break;
			}
			case OGRwkbGeometryType::wkbMultiLineString: {
				for (const auto& p_lineString : *p_geometry->toMultiLineString()) {
					buffOuterPolygons->addGeometry(p_lineString->Buffer(buffOuter));
					if (buffInner != 0) {
						buffInnerPolygons->addGeometry(p_lineString->Buffer(buffInner));
					}
				}
				break;
			}
			default: {
				throw std::runtime_error("geometry type must be LineString or MultiLineString");
			}
		}
	}	

	//step 3: generate the polygon mask and free no longer used memory
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
	
	//step 4: create new GDAL dataset to rasterize as access mask
	GDALDataset *p_accessPolygonDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create(
		"",
		0,
		0,
		0,
		GDT_Unknown,
		nullptr
	);

	OGRLayer *p_layer = p_accessPolygonDataset->CreateLayer(
		"access", 
		p_vector->getLayer(layerName.c_str())->GetSpatialRef(), 
		wkbPolygon, 
		nullptr
	);
	OGRFieldDefn field("index", OFTInteger);
	p_layer->CreateField(&field);
	
	//step 6: add the access polygon to the new layer
	OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
	p_feature->SetField("index", 0);
	p_feature->SetGeometry(p_polygonMask);
	p_layer->CreateFeature(p_feature); //error handling here???
	OGRFeature::DestroyFeature(p_feature);

	//step 7: create in-memory dataset
	GDALDataset *p_accessRasterDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create(
		"",
		p_raster->getWidth(),
		p_raster->getHeight(),
		0,
		GDT_Byte,
		nullptr
	);

	//allocate new raster layer
	void *datapointer = CPLMalloc(p_raster->getWidth() * p_raster->getHeight()* sizeof(uint8_t));

	//add band to new in-memory raster dataset
	char **papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", std::to_string((size_t)datapointer).c_str());

	//step 8: set and fill parameters of new in-memory dataset
	p_accessRasterDataset->AddBand(GDT_Byte, papszOptions);
	p_accessRasterDataset->SetGeoTransform(p_raster->getGeotransform());
	p_accessRasterDataset->SetProjection(p_raster->getDataset()->GetProjectionRef());
	GDALRasterBand *p_band = p_accessRasterDataset->GetRasterBand(1);
	p_band->SetDescription("access_mask");
	p_band->Fill(1);

	//step 9: generate options list for rasterization	
	char **argv;

	//specify invert rasterization and ALL_TOUCHED true
	//this ensures pixels whos upper-left corner is outside the
	//accessable area don't accidentally get included.
	//The upper left corner is where the geotransform applies to.
	argv = CSLAddString(argv, "-i");
	argv = CSLAddString(argv, "-at");

	//specify the burn value for the polygon
	argv = CSLAddString(argv, "-burn");
	argv = CSLAddString(argv, std::to_string(0).c_str());

	//specify the layer
	argv = CSLAddString(argv, "-l");
	argv = CSLAddString(argv, "access");

	GDALRasterizeOptions *options = GDALRasterizeOptionsNew(argv, nullptr);

	//step 10: rasterize vector creating in-memory dataset
	GDALRasterize(
		nullptr,
		p_accessRasterDataset,
		p_accessPolygonDataset,
		options,
		nullptr
	);

	//step 11: free dynamically allocated rasterization options
	GDALRasterizeOptionsFree(options);

	//step11: free no longer used polygon dataset
	free(p_accessPolygonDataset);

	//step 12: return access mask data and pointer to dataset to free
	return {p_accessRasterDataset, datapointer};
}
