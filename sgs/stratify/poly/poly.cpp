/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "helper.h"
#include "raster.h"
#include "vector.h"
#include "gdal_utils.h"

/**
 *
 */
GDALRasterWrapper *poly(
	GDALVectorWrapper *p_vector,
	GDALRasterWrapper *p_raster,
	size_t numStrata,
	std::string layerName,
	std::string query,
	std::string filename,
	bool largeRaster,
	std::string tempFolder,
	std::map<std::string, std::string> driverOptions)
{
	std::cout << "HERE 1" << std::endl;
	GDALAllRegister();

	//step 1: get required info from vector and raster objects
	GDALDataset *p_vectorDS = p_vector->getDataset();
	GDALDataset *p_rasterDS = p_raster->getDataset();
	
	int width = p_raster->getWidth();
	int height = p_raster->getHeight();
	double *geotransform = p_raster->getGeotransform();
	std::string projection = std::string(p_rasterDS->GetProjectionRef());
	std::cout << "HERE 2" << std::endl;

	const OGRSpatialReference *p_layerSrs = p_vector->getLayer(layerName)->GetSpatialRef();
	if (!p_layerSrs) {
		throw std::runtime_error("vector layer does not have a projection.");
	}
	if (projection == "") {
		throw std::runtime_error("raster does not have a projection.");
	}
	if (!OGRSpatialReference(projection.c_str()).IsSame(p_layerSrs)) {
		throw std::runtime_error("raster and vector projections don't match.");
	}
	std::cout << "HERE 3" << std::endl;

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";
	
	RasterBandMetaData band;
	setStratBandTypeAndSize(numStrata - 1, &band.type, &band.size);
	p_rasterDS->GetRasterBand(1)->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
	band.name = "strata";
	std::vector<VRTBandDatasetInfo> VRTBandInfo;
	std::cout << "HERE 4" << std::endl;

	//step 2: create dataset
	GDALDataset *p_dataset = nullptr;
	if (isMEMDataset) {
		p_dataset = createVirtualDataset("MEM", width, height, geotransform, projection);
	       	addBandToMEMDataset(p_dataset, band);
	}
	else if (isVRTDataset) {
		p_dataset = createVirtualDataset("VRT", width, height, geotransform, projection);
		createVRTBandDataset(p_dataset, band, tempFolder, layerName + ".tif", VRTBandInfo, driverOptions); 
	}
	else {
		std::string driver;
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();

		if (extension == ".tif") {
			driver = "Gtiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now");
		}

		p_dataset = createDataset(filename, driver, width, height, geotransform, projection, &band, 1, false, driverOptions);
	}
	std::cout << "HERE 5" << std::endl;

	band.p_band->Fill(band.nan);

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

	std::cout << "HERE 6" << std::endl;
	GDALRasterizeOptions *options = GDALRasterizeOptionsNew(argv, nullptr);
	std::cout << "HERE 6.1" << std::endl;
	//step 5: rasterize vector to in-memory dataset
	GDALRasterize(
		nullptr,
		!isVRTDataset ? p_dataset : VRTBandInfo[0].p_dataset,
		p_vectorDS,
		options,
		nullptr
	);
	std::cout << "HERE 7" << std::endl;

	//step 6: free dynamically allocated rasterization options
	GDALRasterizeOptionsFree(options);
	
	std::cout << "HERE 8" << std::endl;
	if (isVRTDataset) {
		GDALClose(VRTBandInfo[0].p_dataset);
		addBandToVRTDataset(p_dataset, band, VRTBandInfo[0]);	
	}
	std::cout << "HERE 9" << std::endl;

	//step 7: create new GDALRasterWrapper using dataset pointer
	//this dynamically allocated object will be cleaned up by python
	return isMEMDataset ?
		new GDALRasterWrapper(p_dataset, {band.p_buffer}) :
		new GDALRasterWrapper(p_dataset);
}

PYBIND11_MODULE(poly, m) {
	m.def("poly_cpp", &poly);
}
