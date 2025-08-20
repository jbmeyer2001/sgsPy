/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"

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
	size_t splitIndex;
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
	std::vector<void *>stratRasterBands;
	size_t stratRasterBandSize = p_raster->getWidth() * p_raster->getHeight() * sizeof(float);
	void *p_stratRaster = CPLMalloc(stratRasterBandSize * (bandCount + (size_t)map));
	for (size_t i = 0; i < bandCount + (size_t)map; i++) {
		stratRasterBands.push_back((void *)((size_t)p_stratRaster + (stratRasterBandSize * i)));
	}

	//step 2 get raster bands from GDALRasterWrapper
	std::vector<void *>rasterBands;
	std::vector<GDALDataType> rasterBandTypes;
	std::vector<std::vector<double>> probabilities;
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	std::vector<double> noDataValues;
	for (auto const& [key, value] : userProbabilites) {
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
		probabilities.push_back(value);	
		newBandNames.push_back("strat_" + bandNames[key]);
	}

	//step 3 set strat multipliers for mapped stratum raster if required
	std::vector<size_t> bandStratMultipliers(probabilities.size(), 1);
	if (map) {
		for (int i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (probabilities[i].size() + 1);
		}

		newBandNames.push_back("strat_map");
	}
	
	//TODO this can all be done in parallel without much use of locks	
	//step 4 iterate through rasters
	std::vector<std::vector<double>> quantileVals;
	for (int i = 0; i < bandCount; i++) {
		//CALL calculateQuantiles
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
	float noDataFloat = std::nan("-1");
	for (size_t j = 0; j < pixelCount; j++) {
		float mappedStrat = 0;
		float strat = 0;

		for (int i = 0; i < bandCount; i++) {
			std::vector<double> quantiles = quantileVals[i];

			double val;
			switch(rasterBandTypes[i]) {
				case GDT_Int8:
					val = static_cast<double>(reinterpret_cast<int8_t *>(rasterBands[i])[j]);
					break;
				case GDT_UInt16:
					val = static_cast<double>(reinterpret_cast<uint16_t *>(rasterBands[i])[j]);
					break;
				case GDT_Int16:
					val = static_cast<double>(reinterpret_cast<int16_t *>(rasterBands[i])[j]);
					break;
				case GDT_UInt32:
					val = static_cast<double>(reinterpret_cast<uint32_t *>(rasterBands[i])[j]);
					break;
				case GDT_Int32:
					val = static_cast<double>(reinterpret_cast<int32_t *>(rasterBands[i])[j]);
					break;
				case GDT_Float32:
					val = static_cast<double>(reinterpret_cast<float *>(rasterBands[i])[j]);
					break;
				case GDT_Float64:
					val = static_cast<double>(reinterpret_cast<double *>(rasterBands[i])[j]);
					break;
				default:
					throw std::runtime_error("GDALDataType not supported.");
			}

			if (std::isnan(val) || val == noDataValues[i]) {
				strat = noDataFloat;
				mappedStrat = noDataFloat;
			}
			else {	
				strat = static_cast<float>(std::distance(
					quantiles.begin(), 
					std::lower_bound(quantiles.begin(), quantiles.end(), val)
				));
			}

			((float *)stratRasterBands[i])[j] = strat;
			
			if (map) {
				//mappedStrat will remain nodata if it is supposed to be nodata
				mappedStrat += strat * bandStratMultipliers[i] * (mappedStrat != noDataFloat);
			}
		}

		if (map) {
			((float *)stratRasterBands[bandCount])[j] = mappedStrat;
		}
	}

	//free allocated band data
	for (size_t i = 0; i < rasterBands.size(); i++) {
		free(rasterBands[i]);
	}

	//step 9 create new GDALRasterWrapper in-memory
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

	//Step 10 write raster if desired
	if (filename != "") {
		stratRaster->write(filename);
	}

	return stratRaster;
}

PYBIND11_MODULE(quantiles, m) {
	m.def("quantiles_cpp", &quantiles);
}
