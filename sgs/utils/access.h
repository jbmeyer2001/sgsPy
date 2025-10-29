/******************************************************************************
 *
 * Project: sgs
 * Purpose: generate an access mask using vector and raster datasets
 * Author: Joseph Meyer
 * Date: October, 2025
 *
 ******************************************************************************/

#include "helper.h"
#include "raster.h"
#include "vector.h"
#include "gdal_utils.h"

/*
 * This function is old and will no longer be used, and can be removed,
 * once sample_strat has been updated to no longer use it.
 */
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
	char **argv = nullptr;

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

/**
 * This Struct controls the creation and storage of access networks
 * for use in sampling functions.
 */
struct Access {
	bool used = false;
	double area = -1;
	GDALDataset *p_dataset = nullptr;
	RasterBandMetaData band;

	/**
	* This constructor is responsible for setting the used, area, p_dataset, and band
	* members of this struct. In the case where an access vector is given, the dataset
	* which contain a raster dataset with a rasterized version of the access, where a
	* pixel is '1' if it falls within accessible area.
	*
	* First, if p_vector is not given, then the 'used' member remains false.
	*
	* If p_vector is given, it is checked to ensure it has the same spatial reference
	* system as the p_raster. Polygons are created such that The polygons contain
	* the accessible area. This is done by buffering the linestrings in the access
	* vector using buff_outer, and removing the buff_inner buffer.
	*
	* A vector is created usign the output polygons from this calculation. Then,
	* GDALRasterize() is called specifying:
	* -burn 1 				(setting the burn value to 1)
	* -l access 				(specifying the layer name in the polygon dataset as 'access')
	* -te {xmin} {ymin} {xmax} {ymax}	(setting the extent to the raster extent)
	* -ts {width} {height}			(setting the dimensions to the raster dimensions)
	* -ot Int8				(setting the output type to int8_t)
	*
	* The resulting raster dataset created by the GDALRasterize() function is then
	* checked by sampling functions to ensure their samples fall within accessible
	* areas. The 'band' parameter's metadata is set according to the output raster
	* dataset.
	*
	* @param GDALVectorWrapper *p_vector
	* @param GDALRasterWrapper *p_raster
	* std::string layerName
	* double buffInner
	* double buffOuter
	* bool largeRaster
	* std::string tempFolder
	* int xBlockSize
	* int yBlockSize
	*/
	Access(GDALVectorWrapper *p_vector,
		GDALRasterWrapper *p_raster,
		std::string layerName, 
		double buffInner, 
		double buffOuter,
		bool largeRaster,
		std::string tempFolder,
		int xBlockSize,
		int yBlockSize) 
	{
		if (!p_vector) {
			return;
		}

		std::string rastProj = p_raster->getDataset()->GetProjectionRef();
		OGRSpatialReference rastSRS;
		rastSRS.importFromWkt(rastProj.c_str());
		OGRLayer *p_inputLayer = p_vector->getLayer(layerName.c_str());
		OGRSpatialReference *p_vectSRS = p_inputLayer->GetSpatialRef();
		if (!rastSRS.IsSame(p_vectSRS)) {
			throw std::runtime_error("access vector and raster do not have the same spatial reference system.");
		}

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
					OGRGeometry *p_outer = p_geometry->Buffer(buffOuter);
					buffOuterPolygons->addGeometry(p_outer);
	
					if (buffInner != 0) {
						OGRGeometry *p_inner = p_geometry->Buffer(buffInner);
						buffInnerPolygons->addGeometry(p_inner);
					}
	
					break;
				}
				case OGRwkbGeometryType::wkbMultiLineString: {
					for (const auto& p_lineString : *p_geometry->toMultiLineString()) {
						OGRGeometry *p_outer = p_lineString->Buffer(buffOuter);
						buffOuterPolygons->addGeometry(p_outer);
	
						if (buffInner != 0) {
							OGRGeometry *p_inner = p_lineString->Buffer(buffInner);
							buffInnerPolygons->addGeometry(p_inner);
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
	
		this->area = 0;
		switch(wkbFlatten(p_polygonMask->getGeometryType())) {
			case OGRwkbGeometryType::wkbPolygon: {
				OGRCurvePolygon *p_cp = p_polygonMask->toPolygon()->toUpperClass();
				this->area += p_cp->get_Area();
				break;
			}
			case OGRwkbGeometryType::wkbMultiPolygon: {
				for (const auto&p_polygon : *p_polygonMask->toMultiPolygon()) {
					OGRCurvePolygon *p_cp = p_polygon->toUpperClass();
					this->area += p_cp->get_Area();
				}
				break;
			}
			default:
				throw std::runtime_error("p_polygonMask is not a polygon or a multipolygon!");
		}
	
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
	
		std::filesystem::path path = tempFolder;
		path = path / "access.tif";

		//step 9: generate options list for rasterization	
		char **argv = nullptr;
	
		//specify the burn value for the polygon
		argv = CSLAddString(argv, "-burn");
		argv = CSLAddString(argv, std::to_string(1).c_str());
	
		//specify the layer
		argv = CSLAddString(argv, "-l");
		argv = CSLAddString(argv, "access");

		//specify extent
		argv = CSLAddString(argv, "-te");
		argv = CSLAddString(argv, std::to_string(p_raster->getXMin()).c_str());
		argv = CSLAddString(argv, std::to_string(p_raster->getYMin()).c_str());
		argv = CSLAddString(argv, std::to_string(p_raster->getXMax()).c_str());
		argv = CSLAddString(argv, std::to_string(p_raster->getYMax()).c_str());
	
		//specify dimensions
		argv = CSLAddString(argv, "-ts");
		argv = CSLAddString(argv, std::to_string(p_raster->getWidth()).c_str());
		argv = CSLAddString(argv, std::to_string(p_raster->getHeight()).c_str());

		//specify output type
		argv = CSLAddString(argv, "-ot");
		argv = CSLAddString(argv, "Int8");

		GDALRasterizeOptions *options = GDALRasterizeOptionsNew(argv, nullptr);

		//step 10: rasterize vector creating in-memory dataset
		this->p_dataset = GDALDataset::FromHandle(GDALRasterize(
			path.string().c_str(),
			nullptr,
			p_accessPolygonDataset,
			options,
			nullptr
		));

		this->band.p_band = this->p_dataset->GetRasterBand(1);
		this->band.type = GDT_Int8;
		this->band.size = sizeof(int8_t);
		this->band.p_band->GetBlockSize(&this->band.xBlockSize, &this->band.yBlockSize);

		//step 11: free dynamically allocated data
		GDALRasterizeOptionsFree(options);
		free(p_accessPolygonDataset);
	
		this->used = true;
	}
};


