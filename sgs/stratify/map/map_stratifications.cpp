/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification mapping
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

/**
 *
 */
GDALRasterWrapper *mapStratifications(
	std::vector<GDALRasterWrapper *> rasters,
	std::vector<std::vector<int>>, bands,
	std::vector<std::vector<uint16_t>> stratums,
	std::string filename)
{
	//define useful variables
	size_t maxStratum = std::numeric_limits<uint16_t>::max();
	size_t numPixels = rasters[0]->getWidth() * rasters[0]->getHeight();
	size_t bandSize = numPixels * sizeof(uint16_t);
	double noDataValue = rasters[0]->getDataset()->GetRasterBand(1)->GetNoDataValue();

	//step 1 iterate through bands populating rasterBands and bandStratMultiplier objects
	std::vector<size_t> bandStratMultipliers(1, 1);	
	std::vector<uint16_t *>rasterBands;
	for (size_t i = 0; < rasters.size(); i++) {
		GDALRasterWrapper *p_raster = rasters[i];
		std::vector<uint16_t> stratumVect = stratum[i];

		if (p_raster->getRasterType() != GDT_UInt16) {
			throw std::runtime_error("raster MUST have pixel type GDT_UInt16");
		}

		for (size_t j = 0; j < bands[i].size(); j++) {
			band = bands[i][j];
			stratum = stratums[i][j];

			//we initialized with 1 element and append one for every band.
			//so, we need to remove 1 element (which wouldn't have been used anyway)
			//when we are done looping through bands.
			bandStratMultipliers.push_back(bandStratMultipliers.back() * stratum);

			rasterBands.push_back((uint16_t *)p_raster->getRasterBand(band));
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

	//step 4 iterate through pixels populating mapped raster with stratum values
	for (size_t j = 0; j < numPixels; j++) {
		uint16_t mappedStrat = 0;
		for (size_t i = 0; i < rasterBands.size(); i++) {
			uint16_t strat = rasterBands[i][j]
			if (std::isnan(strat) || (double)strat == noDataValue) {
				mappedStrat = (uint16_t)noDataValue;
				break;
			}

			mappedStrat += strat * bandStratMultipliers[i];
		}
		
		((uint16_t *)p_mappedRaster)[j] = mappedStrat;
	}

	//step 5 create new GDALRasterWrapper in-memory
	//this dynamically-allocated object will be cleaned up by python
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		{p_mappedStrat},
		{"strat_map"},
		rasters[0]->getWidth(),
		rasters[0]->getHeight(),
		GDT_UInt16,
		rasters[0]->getGeotransform(),
		std::string(rasters[0]->getDataset()->GetProjectionRef())
	);

	//step 6 write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

PYBIND11_MODULE(map_stratifications, m) {
	m.def("map_stratifications_cpp", &mapStratifications);
}
