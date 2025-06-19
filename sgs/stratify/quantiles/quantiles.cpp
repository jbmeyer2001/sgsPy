/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"

/**
 *
 */
template <typename T, typename U>
GDALRasterWrapper *quantiles(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userProbabilites,
	bool map,
	std::string filename) 
{
	size_t maxStratum = std::numeric_limits<uint16_t>::max();
	int bandCount = userProbabilites.size();

	//step 1 allocate new rasters
	std::vector<void *>stratRasterBands;
	size_t stratRasterBandSize = p_raster->getWidth() * p_raster->getHeight() * sizeof(uint16_t);
	void *p_stratRaster = CPLMalloc(stratRasterBandSize * (bandCount + (size_t)map);
	for (size_t i = 0; i < bandCount + (size_t)map; i++) {
		stratRasterBands.push_back((void *)p_stratRaster + (stratRasterBandSize * i));
	}

	//step 2 get raster bands from GDALRasterWrapper
	std::vector<T *>rasterBands;
	std::vector<std::vector<double>> probabilites
	for (auto const& [key, value] : userProbabilites) {
		rasterBands.push_back((T *)p_raster->getRasterBand(key));

		if (value.size() > maxStratum) {
			throw std::runtime_error("too many stratum given, result would overflow");
		}

		probabilities.push_back(value);	
	}

	//step 3 set strat multipliers for mapped stratum raster if required
	std::vector<size_t> bandStratMultipliers(probabilities.size(), 1);
	if (map) {
		for (int i = 1; i < bandCount; i++) {
			bandStratMultipliers[i] = bandStratMultipliers[i - 1] * (probabilities[i] + 1);
		}
	}
	
	//TODO this can all be done in parallel without much use of locks	
	//step 4 iterate through rasters
	std::vector<std::vector<std::tuple<T, U, uint16_t>>> stratVects(bandCount, {});
	double noDataValue = p_Raster->getDataset->GetRasterBand(1)->GetNoDataValue();
	for (int i = 0; i < bandCount; i++) {
		std::vector<std::tuple<T, U, uint16_t>> *stratVect = &stratVects[i];
		for (size_t j = 0; j < p_raster->getWidth() * p_raster->getHeight(); j++) {
			T val = rasterBands[i][j];

			if (std::isnan(val) || (double)val == noDataValue) {
				//step 4.1 write nodata in nodata ares
				((uint16_t *)stratRasterBands[i])[j] = (uint16_t)noDataValue;
				continue;
			}
			else {
				//step 4.2 populate vectors with index/value
				stratVect->push_back({val, j, 0})
			}
		}	

		//step 4.3 sort vectors in ascending order by pixel val
		std::sort(stratVect->begin(), stratVect->end(), [](auto const& t1, auto const& t2){
			return std::get<0>(t1) < std::get<0>(t2)	
		});
		
		//step 4.4 assign stratum values depending on user defined probability breaks
		size_t numDataPixels = stratVect->size();
		std::vector<auto> splittingIterators;		
		splittingIterators.push_back(stratVect->begin());
		for (const double& prob : probabilities[i]) {
			//TODO test for off by 1!!!!!
			auto it = stratVect->begin() + (size_t)((double)numDataPixels * prob) + 1;
			splittingIndexes.push_back(it);
		}
		splittingIterators.push_back(stratVect->end());

		/**
		 * TODO its possible runtime will be faster if I just write the stratum values here
		 * directly. However, my aim is to take advantage of spatial locality -- to re-sort
		 * the pixels by index then assign them sequentially. For large images, I figure 
		 * the advantage of cache hits will yield a higher runtime despite the additional sort
		 * and memory used to store the stratum value.
		 *
		 * However, this is just a guess, and it may depend on the specific image.
		 *
		 * Regardless, this requires some testing in the future.
		 */
		for (uint16_t i = 0; i < splittingIterators.size() - 1; i++) {
			std::for_each(splittingIterators[i], splittingIterators[i+1], [](auto const& t){
				std::get<3>t = i;		
			});
		}

		//step 4.5 sort by index
		std::sort(stratVect->begin(), stratVect->end(), [](auto const& t1, auto const& t2){
			return std::get<1>(t1) < std::get<1>(t2);		
		});
	}

	//TODO this may be parallelizable
	//step 8 iterate through vectors and assign stratum value
	for (size_t j = 0; j < stratVects[0].size() < j++) {
		uint16_t mappedStrat = 0;
		for (int i = 0; i < bandCount; i++) {
			uint16_t strat = std::get<2>(stratVects[i][j]);
			U index = std::fet<1>(stratVects[i][j]);
			stratRasterBands[i][index] = strat;
			if (map) {
				mappedStrat += strat * bandStratMultipliers[i];
			}
		}

		if (map) {
			stratRasterBands[bandCount][j] = mappedStrat;
		}
	}

	//step 9 create new GDALRasterWrapper in-memory
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	for (const std::string& name : bandNames) {
		newBandNames.push_back("strat_" + name);
	}
	if (map) {
		newBandNames.push_back("strat_map");
	}
	//this dynamically-allocated object will be cleaned up by python
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		stratRasterBands,
		newBandNames,
		p_raster->getWidth(),
		p_raster->getHeight(),
		GDT_UInt16,
		p_raster->getGeotransform(),
		std::string(p_raster->getDataset()->GetProjectionRef())
	);

	//Step 10 write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

/**
 *
 */
GDALRasterWrapper *quantilesTypeSpecifier(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userProbabilites,
	bool map,
	std::string filename) 
{
	//TODO add switch with minIndexIntType and more template arguments
	switch(p_raster->getRasterType()) {
		case GDT_Int8:
		return quantiles<int8_t>(p_raster, userProbabilites, map, filename);
		case GDT_UInt16:
		return breaks<uint16_t>(p_raster, userProbabilites, map, filename);
		case GDT_Int16:
		return breaks<int16_t>(p_raster, userProbabilites, map, filename);
		case GDT_UInt32:
		return breaks<uint32_t>(p_raster, userProbabilites, map, filename);
		case GDT_Int32:
		return breaks<int32_t>(p_raster, userProbabilites, map, filename);
		case GDT_Float32:
		return breaks<float>(p_raster, userProbabilites, map, filename);
		case GDT_Float64:
		return breaks<double>(p_raster, userProbabilites, map, filename);
		default:
		throw std::runtime_error("GDATDataType not one of the accepted types.");
	}
}

PYBIND11_MODULE(quantiles, m) {
	m.def("quantiles_cpp", &quantilesTypeSpecifier);
}
