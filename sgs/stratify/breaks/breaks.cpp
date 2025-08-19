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
	std::vector<void *> rasterBands;
	std::vector<GDALDataType> rasterBandTypes;
	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<std::string> newBandNames;
	std::vector<double> noDataValues;
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

		newBandNames.push_back("strat_" + bandNames[key]);
	}

	//step 3: set bandStratMultipliers and check max size if mapped stratification	
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (int i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (bandBreaks[i].size() + 1);
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
			void *p_data = rasterBands[i];
			double val;
			switch (rasterBandTypes[i]) {
				case GDT_Int8:
					val = static_cast<double>(((int8_t *)p_data)[j]);
					break;
				case GDT_UInt16:
					val = static_cast<double>(((uint16_t *)p_data)[j]);
					break;
				case GDT_Int16:
					val = static_cast<double>(((int16_t *)p_data)[j]);
					break;
				case GDT_UInt32:
					val = static_cast<double>(((uint32_t *)p_data)[j]);
					break;
				case GDT_Int32:
					val = static_cast<double>(((int32_t *)p_data)[j]);
					break;
				case GDT_Float32:
					val = static_cast<double>(((float *)p_data)[j]);
					break;
				case GDT_Float64:
					val = ((double *)p_data)[j];
					break;
				default:
					throw std::runtime_error("raster pixel data type not supported.");
			}

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

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaks);
}
