/******************************************************************************
 *
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using breaks
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include <iostream>

#include "raster.h"
#include "helper.h"

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
GDALRasterWrapper *breaks(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	std::string filename)
{
	//find band count and the maximum number of breaks
	int bandCount = breaks.size();

	//determine max strata for different data types
	int maxInt8 = std::numeric_limits<int8_t>::max();
	int maxInt16 = std::numeric_limits<int16_t>::max();

	//band break locations
	std::vector<std::vector<double>> bandBreaks;
	
	//raster bands and each raster band type
	std::vector<void *> rasterBands;
	std::vector<GDALDataType> rasterBandTypes;
	
	//output strat raster bands and type
	std::vector<void *> stratBands;
	std::vector<GDALDataType> stratBandTypes;

	//band output names
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;

	//per-band nodata values
	std::vector<double> noDataValues;

	//step 1: allocate, read, and initialize raster data and breaks information
	for (auto const& [key, val] : breaks) {
		GDALRasterBand *p_band = p_raster->getRasterBand(key);

		//add band type
		GDALDataType type = p_raster->getRasterBandType(key);
		rasterBandTypes.push_back(type);
		
		//allocate and read raster band buffer
		void *p_data = VSIMalloc3(
			p_raster->getHeight(), 
			p_raster->getWidth(), 
			p_raster->getRasterBandTypeSize(key)
		);	
		CPLErr err = p_band->RasterIO(
			GF_Read,		//GDALRWFlag eRWFlag
			0,			//int nXOff
			0,			//int nYOff
			p_raster->getWidth(),	//int nXSize
			p_raster->getHeight(),	//int nYSize
			p_data,			//void *pData
			p_raster->getWidth(),	//int nBufXSize
			p_raster->getHeight(),	//int nBufYSize
			type,			//GDALDataType eBufType
			0,			//int nPixelSpace
			0			//int nLineSpace
		);
		if (err) {
			throw std::runtime_error("error reading raster band from dataset.");
		}
		rasterBands.push_back(p_data);
			
		//add nodata value to vector
		noDataValues.push_back(p_band->GetNoDataValue());

		//and band breaks vector
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);

		//add type to vector and get size depending on the maximum strata value (given by val.size() + 1)
		size_t pixelTypeSize = setStratBandType(val.size() + 1, stratBandTypes);

		//allocate new strat raster using type size
		void *p_strata = VSIMalloc3(
			p_raster->getWidth(),
			p_raster->getHeight(),
			pixelTypeSize
		);
		stratBands.push_back(p_strata);

		newBandNames.push_back("strat_" + bandNames[key]);
	}

	//step 2: set bandStratMultipliers and check max size if mapped stratification	
	std::vector<size_t> bandStratMultipliers(breaks.size(), 1);
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (int i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (bandBreaks[i].size() + 1);
		}

		//add type to vector and get max size for mapped strat
		size_t maxStrata = bandStratMultipliers.back() * (bandBreaks.back().size() + 1);
		size_t pixelTypeSize = setStratBandType(maxStrata, stratBandTypes);
		void *p_strata = VSIMalloc3(
			p_raster->getWidth(),
			p_raster->getHeight(),
			pixelTypeSize
		);
		stratBands.push_back(p_strata);

		newBandNames.push_back("strat_map");
	}	

	//TODO: multithread and consider cache thrashing
	//step 3: iterate through indices and update the stratified raster bands
	float noDataFloat = std::nan("-1");
	for (size_t j = 0; j < p_raster->getWidth() * p_raster->getHeight(); j++) {
		size_t mappedStrat = 0;
		bool mapNan = false;
		for (int i = 0; i < bandCount; i++) {
			double val = getPixelValueDependingOnType(
				rasterBandTypes[i],
				rasterBands[i],
				j
			);
			//std::cout << "val: " << val << std::endl;			
			bool isNan = std::isnan(val) || (double)val == noDataValues[i];
			mapNan |= isNan;

			//calculate strata value if not nan
			size_t strat = 0;
			if (!isNan) {
				std::vector<double> curBandBreaks = bandBreaks[i];
				auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
				strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
			}

			//std::cout << "strat: " << strat << std::endl;
			setStrataPixelDependingOnType(
				stratBandTypes[i],
				stratBands[i],
				j,
				isNan,
				strat
			);

			//adjust mappedStrat as required
			if (map) {
				mappedStrat += strat * bandStratMultipliers[i];
			}
		}
		
		//assign mapped value
		if (map) {
			setStrataPixelDependingOnType(
				stratBandTypes.back(),
				stratBands.back(),
				j,
				mapNan,
				mappedStrat
			);
		}
	}

	//free allocated band data
	for (size_t i = 0; i < rasterBands.size(); i++) {
		free(rasterBands[i]);
	}

	//step 4: create GDALRasterWrapper object from bands
	//this dynamically-allocated object will be cleaned up by python (TODO I hope...)
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		stratBands, 
		newBandNames,
		stratBandTypes,
		p_raster->getWidth(),
		p_raster->getHeight(),
		p_raster->getGeotransform(),
		std::string(p_raster->getDataset()->GetProjectionRef())
	);
		
	//step 5: write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}
	
	return stratRaster;
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaks);
}
