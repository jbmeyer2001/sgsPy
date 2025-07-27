/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"
#include "vector.h"
#include "gdal_utils.h"

GDALRasterWrapper *poly(
	GDALVectorWrapper *p_vector,
	GDALRasterWrapper *p_raster,
	std::string query,
	std::string filename)
{
	//step 1: get required info from vector and raster objects
	GDALDataset *p_vectorDS = p_vector->getDataset();
	GDALDataset *p_rasterDS = p_raster->getDataset();
	double noDataValue = p_rasterDS->GetRasterBand(1)->GetNoDataValue();
	
	//step 2: create in-memory dataset
	GDALAllRegister();
	GDALDataset *p_dataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create(
		"",
		p_raster->getWidth(),
		p_raster->getHeight(),
		0,
		GDT_Float32,
		nullptr
	);

	//allocate new raster layer
	void *datapointer = CPLMalloc(p_raster->getWidth() * p_raster->getHeight() * sizeof(float));
	
	//add band to new in-memory raster dataset
	char **papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", std::to_string((size_t)datapointer).c_str());
	
	//step 3: set and fill parameters of new in-memory dataset
	p_dataset->AddBand(GDT_Float32, papszOptions);
	p_dataset->SetGeoTransform(p_raster->getGeotransform());
	p_dataset->SetProjection(p_rasterDS->GetProjectionRef());
	GDALRasterBand *p_band = p_dataset->GetRasterBand(1);
	p_band->SetDescription("strata");
	p_band->SetNoDataValue(noDataValue);
	p_band->Fill(noDataValue);

	//step 4: generate options list for GDALRasterize()
	char ** argv = nullptr;

	//specify the vector attribute to 'strata', which is created by the sql query
	argv = CSLAddString(argv, "-a");
	argv = CSLAddString(argv, "strata");

	//specify sql query
	argv = CSLAddString(argv, "-sql");
	argv = CSLAddString(argv, query.c_str());

	//specify sql dialect
	argv = CSLAddString(argv, "-dialect");
	argv = CSLAddString(argv, "SQLITE");	

	GDALRasterizeOptions *options = GDALRasterizeOptionsNew(argv, nullptr);

	//step 5: rasterize vector to in-memory dataset
	GDALRasterize(
		nullptr,
		p_dataset,
		p_vectorDS,
		options,
		nullptr
	);

	//step 6: free dynamically allocated rasterization options
	GDALRasterizeOptionsFree(options);
	
	//step 7: create new GDALRasterWrapper using dataset pointer
	//this dynamically allocated object will be cleaned up by python
	GDALRasterWrapper* stratRaster = new GDALRasterWrapper(p_dataset, {datapointer});

	//step 8: write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

PYBIND11_MODULE(poly, m) {
	m.def("poly_cpp", &poly);
}
