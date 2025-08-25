/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification mapping
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"
#include "helper.h"

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
	size_t height = static_cast<size_t>(rasters[0]->getHeight());
	size_t width = static_cast<size_t>(rasters[0]->getWidth());
	for (size_t i = 1; i < rasters.size(); i++) {
		if (static_cast<size_t>(rasters[i]->getHeight()) != height) {
			std::string err = "raster with index " + std::to_string(i) + " has a different height from the raster at index 0.";
			throw std::runtime_error(err);
		}

		if (static_cast<size_t>(rasters[i]->getWidth()) != width) {
			std::string err = "raster with index " + std::to_string(i) + " has a different width from the raster at index 0.";
			throw std::runtime_error(err);
		}
	}
	

	//step 1 iterate through bands populating rasterBands and bandStratMultiplier objects
	std::vector<size_t> bandStratMultipliers(1, 1);	
	std::vector<void *> rasterBands;
	std::vector<GDALDataType> rasterBandTypes;
	std::vector<double> noDataValues;
	for (size_t i = 0; i < rasters.size(); i++) {
		GDALRasterWrapper *p_raster = rasters[i];
		std::vector<float> stratumVect = stratums[i];

		for (size_t j = 0; j < bands[i].size(); j++) {
			int band = bands[i][j];
			float stratum = stratums[i][j];
			GDALRasterBand *p_band = p_raster->getRasterBand(band);
			GDALDataType type = p_raster->getRasterBandType(band);
			size_t size = p_raster->getRasterBandTypeSize(band);

			printTypeWarningsForInt32Conversion(type);
			rasterBandTypes.push_back(type);
	

			//we initialized with 1 element and append one for every band.
			//so, we need to remove 1 element (which wouldn't have been used anyway)
			//when we are done looping through bands.
			bandStratMultipliers.push_back(bandStratMultipliers.back() * stratum);

			void *p_data = VSIMalloc3(
				height, 
				width, 
				size
			);
			CPLErr err = p_band->RasterIO(
				GF_Read,		//GDALRWFlag eRWFlag
				0,			//int nXOff
				0,			//int nYOff
				width,			//int nXSize
				height,			//int nYSize
				p_data,			//void *pData
				width,			//int nBufXSize
				height,			//int nBufYSize
				type, 			//GDALDataType eBufType
				0,			//int nPixelSpace
				0			//int nLineSpace
			);
			if (err) {
				throw std::runtime_error("error reading raster band from dataset.");
			}
			rasterBands.push_back(p_data);
			noDataValues.push_back(p_band->GetNoDataValue());
		}
	}
	size_t maxStrata = bandStratMultipliers.back();
	bandStratMultipliers.pop_back();

	std::vector<GDALDataType> stratBandTypes;
	size_t size = setStratBandType(maxStrata, stratBandTypes);
	GDALDataType type = stratBandTypes.back();

	//step 3 allocate mapped raster
	void *p_mappedRaster = VSIMalloc3(
		height,
		width,
		size
	);

	//step 4 iterate through pixels populating mapped raster with stratum values
	for (size_t j = 0; j < height * width; j++) {
		size_t mappedStrat = 0;
		bool isNan = false;
		for (size_t i = 0; i < rasterBands.size(); i++) {
			int strat = getPixelValueDependingOnType<int>(
				rasterBandTypes[i],
				rasterBands[i],
				j
			);

			isNan = std::isnan(strat) || (double)strat == noDataValues[i];
			if (isNan) {
				break;
			}
			else {
				mappedStrat += strat * bandStratMultipliers[i];
			}
		}
			
		setStrataPixelDependingOnType(
			type,
			p_mappedRaster,
			j,
			isNan,
			mappedStrat
		);
	}

	//free allocated band data
	for (size_t i = 0; i < rasterBands.size(); i++) {
		free(rasterBands[i]);
	}

	//step 5 create new GDALRasterWrapper in-memory
	//this dynamically-allocated object will be cleaned up by python
	GDALRasterWrapper *stratRaster = new GDALRasterWrapper(
		{p_mappedRaster},
		{"strat_map"},
		stratBandTypes,
		width,
		height,
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
