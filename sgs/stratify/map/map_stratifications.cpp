/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification mapping
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"

/**
 * This function maps multiple stratifications.
 *
 * Based on the number of bands passed, and the number of stratifications
 * per band, multipliers are determined such that every unique combination
 * of stratifications within the passed stratification rasters
 * corrosponds to a single stratification in the output mapping.
 *
 * For example, if there were 3 stratification rasters passed, each with 5
 * possible stratum, there would be 5^3 or 125 possible mapped stratifications.
 *
 * A new GDALRasterWrappper object is created and initialized using the
 * mapped raster, and the mapped raster is written to disk if a filename is given.
 *
 * NOTE: stratifications have the data type 'float' in order to accomodate nan values
 *
 * @param std::vector<GDALRasterWrapper *> rasters list
 * @param std::vector<std::vector<int>> bands list
 * @param std::vector<std::vector<float>> number of stratums list
 * @param std::string filename the filename to write to or "" if file shouldn't be written
 * @returns GDALRasterWrapper *pointer to newly created raster mapping
 */
GDALRasterWrapper *mapStratifications(
	std::vector<GDALRasterWrapper *> rasters,
	std::vector<std::vector<int>> bands,
	std::vector<std::vector<float>> stratums,
	std::string filename)
{
	//define useful variables
	size_t maxStratum = std::numeric_limits<float>::max();
	size_t numPixels = rasters[0]->getWidth() * rasters[0]->getHeight();
	size_t bandSize = numPixels * sizeof(float);

	//step 1 iterate through bands populating rasterBands and bandStratMultiplier objects
	std::vector<size_t> bandStratMultipliers(1, 1);	
	std::vector<float *>rasterBands;
	std::vector<double> noDataValues;
	for (size_t i = 0; i < rasters.size(); i++) {
		GDALRasterWrapper *p_raster = rasters[i];
		std::vector<float> stratumVect = stratums[i];

		if (p_raster->getRasterType() != GDT_Float32) {
			throw std::runtime_error("raster MUST have pixel type GDT_Float32");
		}

		for (size_t j = 0; j < bands[i].size(); j++) {
			int band = bands[i][j];
			float stratum = stratums[i][j];

			//we initialized with 1 element and append one for every band.
			//so, we need to remove 1 element (which wouldn't have been used anyway)
			//when we are done looping through bands.
			bandStratMultipliers.push_back(bandStratMultipliers.back() * stratum);

			rasterBands.push_back((float *)p_raster->getRasterBand(band));
			noDataValues.push_back(p_raster->getDataset()->GetRasterBand(band + 1)->GetNoDataValue());
		}
	}

	//step 2 check max stratum, and remove unused bandStratMultiplier
	//bandStratMultipliers.back() contains the largest index possible given all the stratification
	//combinations.
	if (bandStratMultipliers.back() > maxStratum) {
		throw std::runtime_error("the number of possible strata given by the raster bands exceeds the maximum allowed (the maximum index which can occur without integer overflow).");
	}
	bandStratMultipliers.pop_back();

	//step 3 allocate mapped raster
	void * p_mappedRaster = CPLMalloc(bandSize);
	float noDataFloat = std::nan("-1");
	//step 4 iterate through pixels populating mapped raster with stratum values
	for (size_t j = 0; j < numPixels; j++) {
		float mappedStrat = 0;
		for (size_t i = 0; i < rasterBands.size(); i++) {
			float strat = rasterBands[i][j];
			if (std::isnan(strat) || (double)strat == noDataValues[i]) {
				mappedStrat = noDataFloat;
				break;
			}

			mappedStrat += strat * bandStratMultipliers[i];
		}
			
		((float *)p_mappedRaster)[j] = mappedStrat;
	}

	//step 5 create new GDALRasterWrapper in-memory
	//this dynamically-allocated object will be cleaned up by python
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		{p_mappedRaster},
		{"strat_map"},
		rasters[0]->getWidth(),
		rasters[0]->getHeight(),
		GDT_Float32,
		rasters[0]->getGeotransform(),
		std::string(rasters[0]->getDataset()->GetProjectionRef())
	);
	stratRaster->getDataset()->GetRasterBand(1)->SetNoDataValue(std::nan("-1"));

	//step 6 write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

PYBIND11_MODULE(map_stratifications, m) {
	m.def("map_stratifications_cpp", &mapStratifications);
}
