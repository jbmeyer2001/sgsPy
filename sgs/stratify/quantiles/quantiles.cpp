/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#include "raster.h"
#include "helper.h"

#ifdef WIN32 || _WIN32
	#include <daal_win.h>
#else
	#include <daal.h>
#endif

/*
 *
 */
template <typename T>
std::vector<double> calculateQuantiles(void *p_data, size_t pixelCount, std::vector<double> probabilities, double noDataValue) {
	std::vector<T> values(pixelCount);

	//write values which aren't nodata into values matrix
	bool isNan;
	size_t dataPixelIndex = 0;
	for (size_t i = 0; i < pixelCount; i++) {
		T val = reinterpret_cast<T *>(p_data)[i];
		isNan = std::isnan(val) && (double)val != noDataValue;
		values[dataPixelIndex] = val;
		dataPixelIndex += !isNan;
	}
	size_t numDataPixels = dataPixelIndex + 1;
	values.resize(numDataPixels);

	//sort values matrix in ascending order
	std::sort(values.begin(), values.end());

	//add values which occur at quantile probabilities to return vector
	std::vector<double> quantiles(probabilities.size());
	for (size_t i = 0; i < probabilities.size(); i++) {
		double prob = probabilities[i];
		size_t quantileIndex = (size_t)((double)numDataPixels * prob);
		quantiles[i] = (double)values[quantileIndex];
	}

	return quantiles;
}

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
 * the stratifications as raster bands, and the raster is written to disk
 * if a filename is given.
 *
 * NOTE: stratifications have the data type 'float' to accomodate nan values
 *
 * @param GDALRasterWrapper * a pointer to the raster image to stratify
 * @param std::map<int, std::vector<double>> band mapped to user-defined probs
 * @param bool map whether to add a mapped stratification
 * @param std::string filename the filename to write to (if desired)
 */
GDALRasterWrapper *quantiles(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userProbabilites,
	bool map,
	std::string filename) 
{
	int bandCount = userProbabilites.size();
	size_t pixelCount = p_raster->getHeight() * p_raster->getWidth();

	//step 1 allocate new rasters
	std::vector<void *>stratBands;
	std::vector<GDALDataType> stratBandTypes;
	std::vector<void *>rasterBands;
	std::vector<GDALDataType> rasterBandTypes;
	std::vector<std::vector<double>> probabilities;	
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	std::vector<double> noDataValues;
	for (auto const& [key, val] : userProbabilites) {
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
			p_raster->getHeight(),	//int nBUfYSize
			type,			//GDALDataType eBufType
			0,			//int nPixelSpace
			0			//int nLineSpace
		);
		if (err) {
			throw std::runtime_error("error reading raster band from dataset.");
		}

		rasterBands.push_back(p_data);
		noDataValues.push_back(p_band->GetNoDataValue());
		probabilities.push_back(val);	

		size_t pixelTypeSize = setStratBandType(val.size(), stratBandTypes);

		void *p_strata = VSIMalloc3(
			p_raster->getWidth(),
			p_raster->getHeight(),
			pixelTypeSize
		);
		stratBands.push_back(p_strata);

		newBandNames.push_back("strat_" + bandNames[key]);
	}

	//step 3 set strat multipliers for mapped stratum raster if required
	std::vector<size_t> bandStratMultipliers(probabilities.size(), 1);
	if (map) {
		for (int i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (probabilities[i].size() + 1);
		}

		size_t maxStrata = bandStratMultipliers.back() * probabilities.back().size();
		size_t pixelTypeSize = setStratBandType(maxStrata, stratBandTypes);
		void *p_strata = VSIMalloc3(
			p_raster->getWidth(),
			p_raster->getHeight(),
			pixelTypeSize	
		);
		stratBands.push_back(p_strata);

		newBandNames.push_back("strat_map");
	}
	
	//TODO this can all be done in parallel without much use of locks	
	//step 4 iterate through rasters
	std::vector<std::vector<double>> quantileVals;
	for (int i = 0; i < bandCount; i++) {
		std::vector<double> quantiles;
		switch (rasterBandTypes[i]) {
			case GDT_Int8:
				quantiles = calculateQuantiles<int8_t>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			case GDT_UInt16:
				quantiles = calculateQuantiles<uint16_t>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			case GDT_Int16:
				quantiles = calculateQuantiles<int32_t>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			case GDT_UInt32:
				quantiles = calculateQuantiles<uint32_t>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			case GDT_Int32:
				quantiles = calculateQuantiles<int32_t>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			case GDT_Float32:
				quantiles = calculateQuantiles<float>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			case GDT_Float64:
				quantiles = calculateQuantiles<double>(
					rasterBands[i],
					pixelCount,
					probabilities[i],
					noDataValues[i]
				);
				break;
			default:
				throw std::runtime_error("GDALDataType not supported.");
		}
		quantileVals.push_back(quantiles);
	}

	//TODO this may be parallelizable
	//step 8 iterate through vectors and assign stratum value
	for (size_t j = 0; j < pixelCount; j++) {
		size_t mappedStrat = 0;
		size_t strat = 0;
		bool mapNan = false;

		for (int i = 0; i < bandCount; i++) {
			std::vector<double> quantiles = quantileVals[i];

			double val = getPixelValueDependingOnType<double>(
				rasterBandTypes[i],
				rasterBands[i],
				j		
			);
			
			bool isNan = std::isnan(val) || val == noDataValues[i];
			mapNan |= isNan;

			if (!isNan) {
				strat = std::distance(
					quantiles.begin(), 
					std::lower_bound(quantiles.begin(), quantiles.end(), val)
				);
			}

			setStrataPixelDependingOnType(
				stratBandTypes[i],
				stratBands[i],
				j,
				isNan,
				strat
			);	
			
			if (map && !mapNan) {
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

	//step 9 create new GDALRasterWrapper in-memory
	//this dynamically-allocated object will be cleaned up by python
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		stratBands,
		newBandNames,
		stratBandTypes,
		p_raster->getWidth(),
		p_raster->getHeight(),
		p_raster->getGeotransform(),
		std::string(p_raster->getDataset()->GetProjectionRef())
	);

	//Step 10 write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

PYBIND11_MODULE(quantiles, m) {
	m.def("quantiles_cpp", &quantiles);
}
