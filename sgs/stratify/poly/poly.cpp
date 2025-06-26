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

GDALRasterWrapper *poly(
	GDALVectorWrapper *p_vector,
	GDALRasterWrapper *p_raster,
	std::string attribute,
	std::vector<std::vector<std::string>> features,
	std::string filename)
{
	//generate OGREnvelope
	OGREnvelope envelope;
	envelope.MinX = p_raster->getXMin();
	envelope.MaxX = p_raster->getXMax();
	envelope.MinY = p_raster->getYMin();
	envelope.MaxY = p_raster->getYMax();

	//generate GDALRasterizeOptions using given raster and parameters
	GDALRasterizeOptions* options = GDALRasterizeOptionsNew();
	options->anBandList = {1};
	options->dfXRes = p_raster->getGeotransform[1]; //pixel width 
	options->dfYRes = p_raster->getGeotransform[5]; //pixel height
	options->nXSize = p_raster->getWidth();
	options->nYSize = p_raster->getHeight();
	options->eOutputType = GDT_Float32;
	options->sEnvelop = envelope;
	options->oOutputSRS = p_raster->getDataset()->GetProjectionRef();

	//options to add potentially:
	/*
    	std::vector<double> adfBurnValues{};
    	std::string osFormat{};
    	std::vector<std::string> aosLayers{};
    	std::string osBurnAttribute{};
    	CPLStringList aosRasterizeOptions{};
    	CPLStringList aosTO{};
    	CPLStringList aosCreationOptions{};
    	std::vector<double> adfInitVals{};
    	std::string osNoData{};
    	bool bTargetAlignedPixels = false;
    	bool bCreateOutput = false;
    	*/

	//create new in-memory dataset using driver
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
	GDALDataset *p_dataset = p_driver->Create(
		"",
		options->nXSize,
		options->nYSize,
		options->anBandList.size(),
		options->eOutputType,
		nullptr
	);

	//rasterize image with new in-memory dataset as output
	GDALRasterize(
		nullptr,
		p_dataset,
		p_vector->getDataset(),
		options
	);

	//free rasterize options
	GDALRasterizeOptionsFree(options);

	//create new GDALRasterWrapper using dataset pointer
	//this dynamically allocated object will be cleaned up by python
	GDALRasterWrapper* stratRaster = new GDALRasterWrapper(p_dataset);

	//write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}
}

PYBIND11_MODULE(poly, m) {
	m.def("poly_cpp", &poly);
}
