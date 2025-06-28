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
	double xRes = p_raster->getPixelWidth();
	double yRes = p_raster->getPixelHeight();
	double xMin = p_raster->getXMin();
	double xMax = p_raster->getXMax();
	double yMin = p_raster->getYMin();
	double yMax = p_raster->getYMax();

	//step 2: generate options list
	char ** argv = nullptr;

	//set the attribute to 'strata', which is created by the sql query
	argv = CSLAddString(argv, "-a");
	argv = CSLAddString(argv, "strata");

	//specify sql query
	argv = CSLAddString(argv, "-sql");
	argv = CSLAddString(argv, query.c_str());

	//specify sql dialect
	argv = CSLAddString(argv, "-dialect");
	argv = CSLAddString(argv, "SQLITE");

	//specify nodata
	argv = CSLAddString(argv, "-a_nodata");
	argv = CSLAddString(argv, std::to_string(noDataValue).c_str());

	//specify intialization of pixels to nodata value
	argv = CSLAddString(argv, "-init");
	argv = CSLAddString(argv, std::to_string(noDataValue).c_str());

	//specify resolution
	argv = CSLAddString(argv, "-tr");
	argv = CSLAddString(argv, std::to_string(xRes).c_str());
	argv = CSLAddString(argv, std::to_string(yRes).c_str());

	//specify data type
	argv = CSLAddString(argv, "-ot");
	argv = CSLAddString(argv, "Float32");

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

	//step 3: rasterize vector creating in-memory dataset
	GDALAllRegister();
	GDALDataset *p_dataset = (GDALDataset *)GDALRasterize("",
		nullptr,
		p_vectorDS,
		options,
		nullptr
	);

	//step 4: free dynamically allocated rasterization options
	GDALRasterizeOptionsFree(options);

	//step 5: set the geotransform and projection of the new dataset
	p_dataset->SetGeoTransform(p_raster->getGeotransform());
	p_dataset->SetProjection(p_rasterDS->GetProjectionRef());

	//step 6: create new GDALRasterWrapper using dataset pointer
	//this dynamically allocated object will be cleaned up by python
	GDALRasterWrapper* stratRaster = new GDALRasterWrapper(p_dataset);

	//step 7: write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

PYBIND11_MODULE(poly, m) {
	m.def("poly_cpp", &poly);
}
