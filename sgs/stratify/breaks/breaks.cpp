/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using breaks
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"

/**
 * This function stratifies a given raster using user-defined breaks.
 * The breaks are provided as a vector of doubles mapped to a band index.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map variable to combine the stratification of
 * the multiple raster bands.
 *
 * The required raster bands are first acquired from the datset, and
 * checked to ensure that the number of breaks fits without overflowing.
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
 * NOTE: stratifications have a data type of 'float' to accomodate nan values
 *
 * @param GDALRasterWrapper * a pointer to the raster image to stratify
 * @param std::map<int, std::vector<double>> band mapped to user-defiend breaks
 * @param bool map whether to add a mapped stratification
 * @param std::string filename the filename to write to (if desired)
 */
template <typename T>
GDALRasterWrapper *breaks(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	std::string filename)
{
	//find band count and the maximum number of breaks
	size_t maxBreaks = std::numeric_limits<float>::max();
	int bandCount = breaks.size();

	//step 1: allocate new stratification raster
	std::vector<size_t> bandStratMultipliers(breaks.size(), 1);
	std::vector<std::vector<double>> bandBreaks;
	std::vector<void *> stratRasterBands;
	size_t stratRasterBandSize = p_raster->getWidth() * p_raster->getHeight() * sizeof(float);
	void *p_stratRaster = CPLMalloc(stratRasterBandSize * (bandCount + (size_t)map));
	for (size_t i = 0; i < bandCount + (size_t)map; i++) {
		stratRasterBands.push_back((void *)((size_t)p_stratRaster + (stratRasterBandSize * i)));
	}
	
	//step 2: get the raster bands needed from the datset
	std::vector<T *> rasterBands;
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	std::vector<double> noDataValues;
	for (auto const& [key, val] : breaks) {
		//add requested band to rasterBands
		rasterBands.push_back((T *)p_raster->getRasterBand(key));		
		
		//add nodata value to vector
		noDataValues.push_back(p_raster->getDataset()->GetRasterBand(key + 1)->GetNoDataValue());

		//and band breaks vectori
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);

		//error checking on band count
		if (maxBreaks < val.size() + 1) {
			throw std::runtime_error("number of break indexes exceeds maximum");
		}

		newBandNames.push_back("strat_" + bandNames[key]);
	}

	//step 3: set bandStratMultipliers and check max size if mapped stratification	
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (int i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (bandBreaks[i].size() + 1);
		}

		if (maxBreaks < bandStratMultipliers.back() * bandBreaks.back().size()) {
			throw std::runtime_error("number of break indexes in mapped stratification exceeds maximum.");
		}

		newBandNames.push_back("strat_map");
	}	

	//TODO: multithread and consider cache thrashing
	//step 4: iterate through indices and update the stratified raster bands
	float noDataFloat = std::nan("-1");
	for (size_t j = 0; j < p_raster->getWidth() * p_raster->getHeight(); j++) {
		float mappedStrat = 0;
		bool nan = false;
		for (int i = 0; i < bandCount; i++) {
			T val = rasterBands[i][j];

			if (std::isnan(val) || (double)val == noDataValues[i]) {
				((float *)stratRasterBands[i])[j] = noDataFloat;
				nan = true;
				continue;
			}

			std::vector<double> curBandBreaks = bandBreaks[i];	
			auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
			float strat = (it == curBandBreaks.end()) ? (float)curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
			((float *)(stratRasterBands[i]))[j] = strat;

			if (map) {
				mappedStrat += strat * bandStratMultipliers[i];
			}
		}
		
		if (map) {
			((float *)stratRasterBands[bandCount])[j] = nan ? noDataFloat : mappedStrat;
		}
	}

	//step 5: create GDALRasterWrapper object from bands
	//this dynamically-allocated object will be cleaned up by python
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		stratRasterBands, 
		newBandNames,
		p_raster->getWidth(),
		p_raster->getHeight(),
		GDT_Float32,
		p_raster->getGeotransform(),
		std::string(p_raster->getDataset()->GetProjectionRef())
	);
	GDALDataset *p_dataset = stratRaster->getDataset();
	for (size_t i = 1; i <= stratRasterBands.size(); i++) {
		p_dataset->GetRasterBand(i)->SetNoDataValue(noDataFloat);
	}
	
	//step 6: write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}
	
	return stratRaster;
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
 * A call ismade to breaks() with the necessary data type template
 * argument.
 *
 * @returns GDALRasterWrapper *stratified raster
 */
GDALRasterWrapper *breaksTypeSpecifier(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userDefinedBreaks,
	bool map,
	std::string filename)
{
	switch(p_raster->getRasterType()) {
		case GDT_Int8:
		return breaks<int8_t>(p_raster, userDefinedBreaks, map, filename);
		case GDT_UInt16:
		return breaks<uint16_t>(p_raster, userDefinedBreaks, map, filename);
		case GDT_Int16:
		return breaks<int16_t>(p_raster, userDefinedBreaks, map, filename);
		case GDT_UInt32:
		return breaks<uint32_t>(p_raster, userDefinedBreaks, map, filename);
		case GDT_Int32:
		return breaks<int32_t>(p_raster, userDefinedBreaks, map, filename);
		case GDT_Float32:
		return breaks<float>(p_raster, userDefinedBreaks, map, filename);
		case GDT_Float64:
		return breaks<double>(p_raster, userDefinedBreaks, map, filename);
		default:
		throw std::runtime_error("GDATDataType not one of the accepted types.");
	}
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaksTypeSpecifier);
}
