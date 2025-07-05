/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of stratified sampling
 * Author: Joseph Meyer
 * Date: June, 2025

 *
 ******************************************************************************/

#include <iostream>
#include <random>

#include "access.h"
#include "raster.h"
#include "vector.h"
#include "write.h"

/**
 *
 */
template <typename U>
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_random(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	std::string allocation,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string filename)
{
	//step 1: get raster band
	if (p_raster->getRasterType() != GDT_Float32) {
		throw std::runtime_error("raster band must be of type float32.");
	}
	float *p_strat = p_raster->getBand(0);

	//step 2: create stratum index storing vectors
	//TODO require that users pass us number of stratum!
	std::vector<std::vector<U>> stratumIndexes;
	//stratumIndexes.resize() call resize with the max stratum as

	//step 3: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using mask
		p_mask = maskInfo.second; //pointer to mask
	}

	//stpe 4: iterate through raster band
	double noDataValue = p_dataset->GetDataset()->GetRasterBand(1)->GetNoDataValue();
	U noDataPixelCount = 0;
	for (U i = 0; i < (U)p_raster->getWidth() * (U)p_raster->getHeight(); i++) {
		if (p_access) {
			//step 4.1 if the current pixel is not accessable, mark it as nodata and don't read it
			if (((uint8_t *)p_mask)[i] == 0) {
				noDataPixelCount++;
				continue;
			}
		}

		float val = p_strat[i];
		if (std::isnan(val) || (double)val == noDataValue) {
			//step 4.2: increment noDataPixelCount if encountered noData
			noDataPixelCount++;
		}
		else {
			//step 4.3: add data index to stratum index array
			stratumIndexes[(size_t)val].push_back(i);
		}
	}
	if (p_access) {
		free(p_accessMaskDataset);
	}

	//step 5: using the size of the stratum and allocation method,
	//determine the number of samples to take from each stratum
	

	//make a map and an unordered map in the same way that srs does
	//
	//generate random number generator
	//
	//add sample pixels
	//
	//calculate geotransform of sample pixels
}
	
/**
 *
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_queinnec(
	GDALRasterWrapper,
	size_t numSamples,
	int wrow,
	int wcol,
	std::string allocation,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string fileName)
{

}

/**
 *
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_cpp_access(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	int wrow,
	int wcol,
	std::string allocation,
	std::string method,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string filename)
{
	if (method == "random") {
		return strat_random(
			p_raster,
			numSamples,
			allocation,
			weights,
			mindist,
			p_access,
			layerName,
			buffInner,
			buffOuter,
			filename
		);
	}
	else { //method == "Queinnec"
		return strat_queinnec(
			p_raster,
			numSamples,
			wrow,
			wcol,
			allocation,
			weights,
			mindist,
			p_access,
			layerName,
			buffInner,
			buffOuter,
			filename
		);
	}
}

/**
 *
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_cpp(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	int wrow,
	int wcol,
	std::string allocation,
	std::string method,
	std::vector<double> weights,
	double mindist,
	std::string filename)
{
	if (method == "random") {
		return strat_random(
			p_raster,
			numSamples,
			allocation,
			weights,
			mindist,
			nullptr,
			"",
			0,
			0,
			filename
		);
	}
	else { //method == "Queinnec"
		return strat_queinnec(
			p_raster,
			numSamples,
			wrow,
			wcol,
			allocation,
			weights,
			mindist,
			nullptr,
			"",
			0,
			0,
			filename
		);
	}

}

PYBIND11_MODULE(strat, m) {
	m.def("strat_cpp", &strat_cpp);
	m.def("strat_cpp_access", &strat_cpp_access);
}
