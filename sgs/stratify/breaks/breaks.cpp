/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratificaiton using breaks
 * Author: Joseph Meyer
 * Date: June, 2025

 *
 ******************************************************************************/

#include "raster.h"
#include "write.h"

/**
 *
 */
template <typename T>
std::tuple<GDALRasterWrapper *, std::vector<double>, std::map<double, int>>
breaks(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	bool plot)
{
	size_t maxBreaks = std::numeric_limits<uint16_t>::max();
	int bandCount = breaks.size();

	//step 1: get dataset
	GDALDataset *p_dataset = p_raster->getDataset();

	//step 2: allocate new stratification raster usign uint16_t
	std::vector<uint16_t *>stratRasterBands;
	std::vector<uint16_t> bandStratMultipliers(bands.size(), 1);
	std::vector<std::vector<double>> bandBreaks;
	size_t stratRasterLayerSize = p_raster->getWidth() * p_raster->getHeight() * sizeof(uint16_t);
	uint16_t *p_stratRaster = (uint16_t *)CPLMalloc(
		stratRasterLayerSize * 
		(size_t)(bandCount + (int)map)
	);
	
	//step 3: allocate the raster bands
	std::vector<T *> rasterBands;
	size_t rasterLayerSize = p_raster->getWidth() * p_raster->getHeight() * sizeof(T);
	T *p_raster = (T *)CPLMalloc(
		rasterLayerSize *
		(size_t)bandCount
	);
	
	//step 4: read the raster band
	CPLErr err;
	void *p_stratRasterBuffer = (void *)p_stratRaster;
	void *p_rasterBuffer = (void *)p_raster;
	for (auto const& [key, val] : breaks) {
		err = p_dataset->GetRasterBand(key + 1)->RasterIO(
			GF_Read,			//GDALRWFlag eRWFlag
			0,				//int nXOff
			0,				//int nYOff
			p_dataset->getWidth(),		//int nXSize
			p_dataset->getheight(),		//int nYSize
			p_rasterBuffer,			//void *pData
			p_dataset->getWidth(),		//int nBufXSize
			p_dataset->getHeight(),		//int nBufYSize
			p_dataset->getRasterType(),	//GDALDataType eBufType
			0,				//int nPixelType
			0				//int nLineSpace
		);
		if (err) {
			throw std::runtime_error("error reading raster band from dataset.");
		}
			
		//move the write buffer by the size of a band
		p_rasterBuffer = (void *)((size_t)p_rasterBuffer + rasterLayerSize);

		//populate band vectors with pointers
		p_stratRasterBuffer = (void *)((size_t)p_stratRasterBuffer + stratRasterLayerSize);
		rasterBands.append((T *)p_rasterBuffer);
		stratRasterBands.append((uint16_t *)p_stratRasterBuffer);

		//and band breaks
		bandBreaks.append(val);

		//error checking on band count
		if (maxBreaks < val.size() + 1) {
			throw std::runtime_error("number of break indexes (" + std::to_string(val.size() + 1) + ") exceeds maximum of " + std::to_string(maxBreaks) ".");
		}
	}

	//if map is true add an extra map raster band
	if (map) {
		p_stratRasterBuffer = (void *)((size_t)p_stratRasterBuffer + stratRasterLayerSize);
		stratRasterBands.append((uint16_t *)p_stratRasterBuffer);

		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (int i = 1; i < bandCount; i++) {
			bandStratMultipliers[i] = bandStratMultipliers[i - 1] * bandBreaks
		}
	}
	
	//TODO: multithread and consider cache thrashing
	//step 6: iterate through indices
	for (size_t j = 0; j < p_dataset->getWidth() * p_dataset->getHeight(); j++) {
		uint16_t mappedStrat = 0;
		for (int i = 0; i < bandCount; i++) {
			T val = rasterBands[i][j];
			breaks = bandBreaks[i];
			auto upper std::upper_bound(breaks.begin(), breaks.end(), val);
			uint16_t strat = (upper == breaks.end()) ? breaks.size() - 1 : std::distance(breaks.begin(), upper);
			stratRasterBands[i][j] = strat;

			if (map) {
				mappedStrat += strat * bandStratMultipliers[i];
			}

			if (plot) {
				//DO SOMETHING
			}
		}
		
		if (map) {
			stratRasterBands[bandCount][j] = mappedStrat;
		}
	}

	//Step 8: free no longer required memory
	CPLFree(p_raster);

	//step 8: create GDALRasterWrapper object from bands
	//TODO if write is true...
	//driver
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverBYName("MEM");
	GDALDataset *p_stratDataset = p_driver->Create(
		'', 
		p_dataset->getWidth(), 
		p_dataset->getHeight(),
	      	0, //if write is true this should be band count + (int)map
		GDT_UInt16
	);

	//for bands (if write is not true)
	//AddBand(GDT_Uint16, 
	//papszOptions = "DATAPOINTER= pointer to data
	//use CPLPrintPointer ??? 

	return {p_stratRaster, {}, {}}
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
std::tuple<GDALRasterWrapper *, std::vector<double>, std::map<double, int>>
breaksTypeSpecifier(
	GDALRasterWrapper *p_raster,
	std::vector<std::vector<double>> breaks,
	bool map,
	bool plot)
{
	switch(p_raster->getRasterType()) {
		case GDT_Int8:
		return breaks<int8_t>(p_raster, breaks, map, plot);
		case GDT_UInt16:
		return breaks<uint16_t>(p_raster, breaks, map, plot);
		case GDT_Int16:
		return breaks<int16_t>(p_raster, breaks, map, plot);
		case GDT_UInt32:
		return breaks<uint32_t>(p_raster, breaks, map, plot);
		case GDT_Int32:
		return breaks<int32_t>(p_raster, breaks, map, plot);
		case GDT_Float32:
		return breaks<float>(p_raster, breaks, map, plot);
		case GDT_Float64:
		return breaks<double>(p_raster, breaks, map, plot);
		default:
		throw std::runtime_error("GDATDataType not one of the accepted types.");
	}
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaksTypeSpecifier);
}
