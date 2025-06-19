/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using breaks
 * Author: Joseph Meyer
 * Date: June, 2025

 *
 ******************************************************************************/

#include <iostream> //TODO remove

#include "raster.h"
#include "write.h"

/**
 * this function stratifies a given raster using user-defined breaks.
 * The breaks are provided as a vector of doubles corresponding to a band index.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map function, to combine the stratification of
 * the multiple raster bands.
 *
 * The required raster bands are first acquired from the datset, and
 * checked to ensure that the number of breaks fits without overflowing.
 * plotInfo, mins, and bucketSizes are updated if plot has been defined.
 *
 * An additional band is added to the output raster if map is specified,
 * and multipliers are determined, so that every unique combination
 * of stratifications of the normal raster bands corresponds to a 
 * single stratification of the mapped raster band. For example,
 * if there were 3 normal raster bands each with 5 possible stratifications,
 * there would be 5^3 or 125 possible mapped stratificaitons.
 *
 * The raster bands are then iterated thorough and stratifications are determined.
 * a new dataset object is created using the stratifcation raster, and it is
 * written to disk if a filename is given.
 *
 * @param GDALRasterWrapper * a pointer to the raster image were stratifying
 * @param std::map<int, std::vector<double>> band and user-defiend breaks mapping
 * @param bool map whether to add a mapped stratification
 * @param bool plot whether to accumulate and return plotting information
 * @param std::string filename the filename to write to (if desired)
 */
template <typename T>
std::pair<GDALRasterWrapper *, std::vector<std::vector<size_t>>>
breaks(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	bool plot,
	std::string filename)
{
	//define plotting data structures
	std::vector<std::vector<size_t>> plotInfo;
	std::vector<double> mins;
	std::vector<double> bucketSizes;

	//find band count and the maximum number of breaks
	size_t maxBreaks = std::numeric_limits<uint16_t>::max();
	int bandCount = breaks.size();

	//step 1: get dataset
	GDALDataset *p_dataset = p_raster->getDataset();

	//step 2: allocate new stratification raster
	std::vector<size_t> bandStratMultipliers(breaks.size(), 1);
	std::vector<std::vector<double>> bandBreaks;
	std::vector<uint16_t *> stratRasterBands;
	size_t stratRasterLayerSize = p_raster->getWidth() * p_raster->getHeight() * sizeof(float);
	for (int i = 0; i < bandCount + (int)map; i++) {
		stratRasterBands.push_back((uint16_t *)CPLMalloc(stratRasterLayerSize));
	}
	
	//step 3: get the raster bands needed from the datset
	CPLErr err;
	std::vector<void *> rasterBands;
	for (auto const& [key, val] : breaks) {
		//add requested band to rasterBands
		rasterBands.push_back((void *)p_raster->getRasterBand(key));		

		//and band breaks vector
		bandBreaks.push_back(val);

		//error checking on band count
		if (maxBreaks < val.size() + 1) {
			throw std::runtime_error("number of break indexes exceeds maximum");
			//throw std::runtime_error("number of break indexes (" + std::to_string(val.size() + 1) + ") exceeds maximum of " + std::to_string(maxBreaks) ".");
		}

		//add min and bucket size for histogram plot
		if (plot) {
			double min = p_raster->getMinPixelVal(key);
			double max = p_raster->getMinPixelVal(key);
			std::vector<size_t> buckets(30, 0);
			mins.push_back(min);
			bucketSizes.push_back((max - min) / 30); //30 buckets
			plotInfo.push_back(buckets);

		}
	}

	//step 4: set bandStratMultipliers and check max size if mapped stratification	
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (int i = 1; i < bandCount; i++) {
			bandStratMultipliers[i] = bandStratMultipliers[i - 1] * bandBreaks[i].size();
		}

		if (maxBreaks < bandStratMultipliers[bandCount - 1] * bandBreaks[bandCount - 1].size()) {
			throw std::runtime_error("number of break indexes in mapped stratification exceeds maximum");
		}
	}
	
	//TODO: multithread and consider cache thrashing
	//step 5: iterate through indices and update the stratified raster bands
	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
	for (size_t j = 0; j < p_raster->getWidth() * p_raster->getHeight(); j++) {
		uint16_t mappedStrat = 0;
		bool nan = false;
		for (int i = 0; i < bandCount; i++) {
			T val = ((T *)rasterBands[i])[j];
			if (std::isnan(val) || (double)val == noDataValue) {
				stratRasterBands[i][j] = (uint16_t)noDataValue;
				nan = true;
				continue;
			}

			std::vector<double> curBandBreaks = bandBreaks[i];	
			auto upper = std::upper_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
			uint16_t strat = (upper == curBandBreaks.end()) ? (uint16_t)breaks.size() - 1 : std::distance(curBandBreaks.begin(), upper);
			stratRasterBands[i][j] = strat;

			if (map) {
				mappedStrat += strat * bandStratMultipliers[i];
			}

			if (plot) {
				//calculate bucket and increment corresponding bucket info
				size_t bucket = (size_t)((val - mins[i]) / bucketSizes[i]);
				plotInfo[i][bucket]++;
			}
		}
		
		if (map) {
			stratRasterBands[bandCount][j] = nan ? (uint16_t)noDataValue : mappedStrat;
		}
	}

	//step 6: create GDALRasterWrapper object from bands
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	for (int i = 0; i < bandNames.size(); i++) {
		newBandNames.push_back("strat_" + bandNames[i]);
	}
	if (map) {
		newBandNames.push_back("strat_map");
	}
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		rasterBands, 
		newBandNames,
		p_raster->getWidth(),
		p_raster->getHeight(),
		GDT_UInt16, 
		p_raster->getGeotransform(),
		std::string(p_dataset->GetProjectionRef())
	);
	
	//step 7: write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}
	
	return {stratRaster, plotInfo};
}

/**
 * Having template types which rely on dynamic information (such as pixel
 * type of the raster) require an unfortunate amount of boilerplate code.
 *
 * This is an attempt to condense as much of the annoying boilerplate into
 * a single place.
 *
 * This function uses type information of the raster pixel type.
 *
 * A call ismade to breads() with the necessary data type template
 * argument.
 *
 * @returns std::tuple<GDALRasterWrapper *, std::vector<double>, std::map<double, int>>
 * 		stratified raster, and plotting information
 */
std::pair<GDALRasterWrapper *, std::vector<std::vector<size_t>>>
breaksTypeSpecifier(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userDefinedBreaks,
	bool map,
	bool plot,
	std::string filename)
{
	switch(p_raster->getRasterType()) {
		case GDT_Int8:
		return breaks<int8_t>(p_raster, userDefinedBreaks, map, plot, filename);
		case GDT_UInt16:
		return breaks<uint16_t>(p_raster, userDefinedBreaks, map, plot, filename);
		case GDT_Int16:
		return breaks<int16_t>(p_raster, userDefinedBreaks, map, plot, filename);
		case GDT_UInt32:
		return breaks<uint32_t>(p_raster, userDefinedBreaks, map, plot, filename);
		case GDT_Int32:
		return breaks<int32_t>(p_raster, userDefinedBreaks, map, plot, filename);
		case GDT_Float32:
		return breaks<float>(p_raster, userDefinedBreaks, map, plot, filename);
		case GDT_Float64:
		return breaks<double>(p_raster, userDefinedBreaks, map, plot, filename);
		default:
		throw std::runtime_error("GDATDataType not one of the accepted types.");
	}
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaksTypeSpecifier);
}
