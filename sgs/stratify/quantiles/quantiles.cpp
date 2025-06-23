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
 * This function stratifies a given raster using user-defined probabilities.
 * The probabilities (quantiles) are provided as a vector of doubles mapped
 * to a band index.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map variable to combine the stratification of
 * the given raster bands.
 *
 * The raster bands are initially iterated through, and vectors are 
 * creatied containing their value, index in the original raster, and 
 * resulting stratum (not set). The vectors are sorted by value, and
 * their stratum is set according to which quantile they are in.
 *
 * The vectors are then re-sorted by index -- in order to take advantage
 * of the CPU cache and spatial locality for large images -- before being
 * iterated through (now sequentially by index) and having their stratum
 * written to the output raster. 
 *
 * An additional band is added to the output raster if map is specified,
 * and multipliers are determined, so that every unique combination of 
 * stratifications of the normal raster bands corresponds to a single
 * stratification of the mapped raster band. For example, if there were
 * 3 normal raster bands each with 5 possible stratifications, there
 * would be 5^3 or 125 possible mapped stratifications.
 *
 * A new GDALRasterWrapper object is created and initialized using
 * the strtified as raster bands, and the raster is written to disk
 * if a filename is given.
 *
 * @param GDALRasterWrapper * a pointer to the raster image to stratify
 * @param std::map<int, std::vector<double>> band mapped to user-defined probs
 * @param bool map whether to add a mapped stratification
 * @param std::string filename the filename to write to (if desired)
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
	void *p_stratRaster = CPLMalloc(stratRasterBandSize * (bandCount + (size_t)map));
	for (size_t i = 0; i < bandCount + (size_t)map; i++) {
		stratRasterBands.push_back((void *)((size_t)p_stratRaster + (stratRasterBandSize * i)));
	}

	//step 2 get raster bands from GDALRasterWrapper
	std::vector<T *>rasterBands;
	std::vector<std::vector<double>> probabilities;
	std::vector<std::vector<std::tuple<T, U, uint16_t>>> stratVects;
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	for (auto const& [key, value] : userProbabilites) {
		rasterBands.push_back((T *)p_raster->getRasterBand(key));

		if (value.size() > maxStratum) {
			throw std::runtime_error("too many stratum given, result would overflow");
		}

		probabilities.push_back(value);

		std::vector<std::tuple<T, U, uint16_t>> stratVect;
		stratVects.push_back(stratVect);

		newBandNames.push_back("strat_" + bandNames[key]);
	}

	//step 3 set strat multipliers for mapped stratum raster if required
	std::vector<size_t> bandStratMultipliers(probabilities.size(), 1);
	if (map) {
		for (int i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (probabilities[i].size() + 1);
		}
		
		if (maxBreaks < bandStratMultipliers[bandCount - 1] * (probabilities[bandCount - 1].size() + 1)) {
			throw std::runtime_error("number of stratum indexes in mapped stratification exceeds maximum.");
		}

		newBandNames.push_back("strat_map");
	}
	
	//TODO this can all be done in parallel without much use of locks	
	//step 4 iterate through rasters
	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
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
				stratVect->push_back({val, j, 0});
			}
		}	

		//step 4.3 sort vectors in ascending order by pixel val
		std::sort(stratVect->begin(), stratVect->end(), [](auto const& t1, auto const& t2){
			return std::get<0>(t1) < std::get<0>(t2);	
		});
		
		//step 4.4 assign stratum values depending on user defined probability breaks
		size_t numDataPixels = stratVect->size();
		std::vector<size_t> stratumSplittingIndexes;		
		stratumSplittingIndexes.push_back(0);
		for (const double& prob : probabilities[i]) {
			//TODO test for off by 1!!!!!
			stratumSplittingIndexes.push_back((size_t)((double)numDataPixels * prob) + 1);
		}
		stratumSplittingIndexes.push_back(numDataPixels);

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
		for (uint16_t stratum = 0; stratum < stratumSplittingIndexes.size() - 1; stratum++) {
			for (size_t index = stratumSplittingIndexes[stratum]; index < stratumSplittingIndexes[stratum + 1]; index++) {
				std::get<2>(stratVect->at(index)) = stratum;
			}
		}

		//step 4.5 sort by index
		std::sort(stratVect->begin(), stratVect->end(), [](auto const& t1, auto const& t2){
			return std::get<1>(t1) < std::get<1>(t2);		
		});
	}

	//TODO this may be parallelizable
	//step 8 iterate through vectors and assign stratum value
	for (size_t j = 0; j < stratVects[0].size(); j++) {
		uint16_t mappedStrat = 0;

 		//I'm assuming if there are nodata pixels, they're all in the same indexes.
		//if that's not true there could be errors here.
		U index = std::get<1>(stratVects[0][j]);
		for (int i = 0; i < bandCount; i++) {
			uint16_t strat = std::get<2>(stratVects[i][j]);
			((uint16_t *)stratRasterBands[i])[index] = strat;
			if (map) {
				mappedStrat += strat * bandStratMultipliers[i];
			}
		}

		if (map) {
			((uint16_t *)stratRasterBands[bandCount])[index] = mappedStrat;
		}
	}

	//step 9 create new GDALRasterWrapper in-memory
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
 * Having template types which rely on dynamic information (such as
 * the pixel type of the added raster, or the number of pixels in
 * the raster) require an unfortunate amount of boilerplate code.
 *
 * This is an attempt to condense as much of the annoyting boilerplate
 * into a single place.
 *
 * This function uses type information of the raster pixel type,
 * as well as the minimally sized unsigned int type which can represent
 * all necessary indices.
 *
 * A call is made to quantiles() with the necessary data type template
 * arguments depending on raster parameters.
 *
 * @returns GDALRasterWraper * newly generated raster image of stratum
 */
GDALRasterWrapper *quantilesTypeSpecifier(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userProbabilites,
	bool map,
	std::string filename) 
{
	//TODO add switch with minIndexIntType and more template arguments
	std::string minIndexIntType = p_raster->getMinIndexIntType(true); //singleBand = true
	switch(p_raster->getRasterType()) {
		case GDT_Int8:
		if (minIndexIntType == "unsigned_short") { return quantiles<int8_t, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<int8_t, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<int8_t, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<int8_t, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		case GDT_UInt16:
		if (minIndexIntType == "unsigned_short") { return quantiles<uint16_t, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<uint16_t, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<uint16_t, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<uint16_t, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		case GDT_Int16:
		if (minIndexIntType == "unsigned_short") { return quantiles<int16_t, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<int16_t, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<int16_t, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<int16_t, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		case GDT_UInt32:
		if (minIndexIntType == "unsigned_short") { return quantiles<uint32_t, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<uint32_t, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<uint32_t, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<uint32_t, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		case GDT_Int32:
		if (minIndexIntType == "unsigned_short") { return quantiles<int32_t, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<int32_t, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<int32_t, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<int32_t, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		case GDT_Float32:
		if (minIndexIntType == "unsigned_short") { return quantiles<float, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<float, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<float, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<float, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		case GDT_Float64:
		if (minIndexIntType == "unsigned_short") { return quantiles<double, unsigned short>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned") { return quantiles<double, unsigned>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long") { return quantiles<double, unsigned long>(p_raster, userProbabilites, map, filename); }
		if (minIndexIntType == "unsigned_long_long") { return quantiles<double, unsigned long long>(p_raster, userProbabilites, map, filename); }
		break;
		default:
		throw std::runtime_error("GDATDataType not one of the accepted types.");
	}
	throw std::runtime_error("type " + minIndexIntType + " not a valid type.");
}

PYBIND11_MODULE(quantiles, m) {
	m.def("quantiles_cpp", &quantilesTypeSpecifier);
}
