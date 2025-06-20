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
	size_t maxStratum = std::numeric_limits<uint16_t>::max();
	size_t numPixels = rasters[0]->getWidth() * rasters[0]->getHeight();
	size_t bandSize = numPixels * sizeof(uint16_t);
	double noDataValue = rasters[0]->getDataset()->GetRasterBand(1)->GetNoDataValue();

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

			//TODO RETHINK BAND STRAT MULTIPLIERS CALCULATIONS
			//size_t numBands = bandStratMultipliers.size();
			//else {
			//	bandStratMultipliers.push_back(bandStratMultipliers[numBands - 1] * stratum);
			//}

			rasterBands.push_back((uint16_t *)p_raster->getRasterBand(band));
		}
	}

	void * p_mappedRaster = CPLMalloc(bandSize * rasterBands.size());
	std::vector<void *> mappedRasterBands;
	for (size_t i = 0; i < rasterBands.size(); i++) {
		mappedRasterBands.append((void *)((size_t)p_mappedRaster + layerSize * i));
	}	

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
		
		((uint16_t *)mappedRasterBands[i])[j] = mappedStrat;
	}
}

PYBIND11_MODULE(map_stratifications, m) {
	m.def("map_stratifications_cpp", &mapStratifications);
}
