/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of PCA
 * Author: Joseph Meyer
 * Date: January, 2026
 *
 ******************************************************************************/

/**
 * @defgroup dist dist
 * @ingroup utils
 */

#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

namespace sgs {
namespace dist {

/**
 * @ingroup dist
 * This is a helper function for finding the minimum and maximum values
 * within a given raster band.
 *
 * The min and max variables are passed as references, and the min value
 * is initialized to the largest value the data type can have, and the max
 * value is initialized to the smallest value the data type can have.
*
 * @param RasterBandMetaData& band
 * @param int width
 * @param int height
 * @param T& min
 * @param T& max
 */
template <typename T>
void 
findMinMax(
	helper::RasterBandMetaData& band, 
	int width, 
	int height,
	T& min,
	T& max) 
{
	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (height + band.yBlockSize - 1) / band.yBlockSize;

	min = std::numeric_limits<T>::max();
	max = std::numeric_limits<T>::min();
	T nan = static_cast<T>(band.nan);
	void *p_data = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);

	//calculate raster band minimum and maximum values
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_mutex->lock();
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			band.p_mutex->unlock();
			
			//read the band into memory
			rasterBandIO(band, p_data, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, true);

			for (int y = 0; y < yValid; y++) {
				int index = y * width;
				for (int x = 0; x < xValid; x++) {
					T val = reinterpret_cast<T *>(p_data)[index];
					if (val != nan && !std::isnan(val)) {
						min = std::min(min, val);
						max = std::max(max, val);
					}
					index++;
				}
			}	
		}
	}
	
	VSIFree(p_data);
}

/**
 * @ingroup dist
 * This function calculates the bin values for the distribution.
 *
 * It takes miniumum and maximum values as parameters. A vector of double
 * values is created to represent the bins, in addition to a vector of 
 * the data type of the raster pixels. The vector of doubles is the one
 * returned to the user, the vector of the raster pixel type is used
 * to calculate which bin a particular pixel goes into while iterating
 * through the raster. 
 *
 * The reason why there is an additional vector with the rasters data type 
 * is that it will reduce the potentialy overhead of requiring the program to cast 
 * the pixel type to double to calculate the bin for every pixel. However, it also
 * means that for integer values a simple static cast will not work, because it will
 * always truncate the value. If the actual bin stops at 10.1, but the vin says 10 instead
 * of 11, then pixels with value 10 will be placed into the wrong bin. The bin value
 * represents the START value of the bin. The final parameter of the double bin 
 * is the 'cap' and represents the end of the last bini, as is customary. 
 * The raster pixel type bin vector does not have this cap value, as it is not necessary.
 *
 * @param T min
 * @param T max
 * @param int nBins
 * @param GDALDataType type
 * @param std::vector<double>& bins
 * @param std::vector<T>& bins
 */
template <typename T>
void
setBins(
	T min,
	T max,
	int nBins,
	GDALDataType type,
	std::vector<double>& dbins,
	std::vector<T>& tbins)
{
	dbins.resize(nBins + 1);
	tbins.resize(nBins);

	//determine bin values as doubles to display to the user.
	double step = ((double)max - (double)min) / ((double)nBins);
	double cur = (double)min;
	for (int i = 0; i <= nBins; i++) {
		dbins[i] = cur;
		cur += step;
	}
	

	//set bin values as the data type used so less data type conversion is necessary in later code
	cur = (double)min;
	if (type == GDT_Float32 || type == GDT_Float64) {
		for (int i = 0; i < nBins; i++) {
			tbins[i] = static_cast<T>(dbins[i]);
		}
	}
	else {
		tbins[0] = min;
		for (int i = 1; i < nBins; i++) {
			cur += step;
			tbins[i] = std::ceil(cur);
		}
	}
}

/**
 * @ingroup dist
 * This function calculates the distribution of pixels across the entire population
 * (the whole raster band).
 *
 * The chunk of the raster assigned to this call of the function is iterated through
 * in blocks. All of the pixels in each block are then iterated through, the nan values
 * are ignored.
 *
 * The bin values are passed as a parameter, and the value in the vector represents the
 * start of a particular bin. This means that the first index in the bins vector which
 * is less than or equal to a particular pixel, is that pixles bin. For each pixel,
 * the iteration starts at the largest bin and iterates down until it reaches a bin
 * value less than or equal to a particular pixel.
 *
 * @param RasterBandMetaData& band
 * @param int width
 * @param int height
 * @param int nBins
 * @param std::vector<T>& binVals
 * @param std::vector<int64_t>& counts
 */
template <typename T>
void
populationDistribution(
	helper::RasterBandMetaData& band,
	int width,
	int height,
	int nBins,
	std::vector<T>& binVals,
	std::vector<int64_t>& counts)
{
	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (height + band.yBlockSize - 1) / band.yBlockSize;

	T nan = static_cast<T>(band.nan);
	void *p_data = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			
			int xValid, yValid;
			band.p_mutex->lock();
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			band.p_mutex->unlock();

			//read the band into memory
			rasterBandIO(band, p_data, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, true);

			for (int y = 0; y < yValid; y++) {
				int index = y * width;
				for (int x = 0; x < xValid; x++) {
					T val = reinterpret_cast<T *>(p_data)[index];

					if (val != nan && !std::isnan(val)) {
						for (int i = nBins - 1; i >= 0; i--) {
							if (binVals[i] <= val) {
								counts[i]++;	
								break;
							}
						}
					}
					index++;
				}
			}
		}
	}
	
	VSIFree(p_data);
}

/**
 * @ingroup dist
 *
 * This function calculates the probability distribution of a sample network, if a vector layer
 * of samples was provided. It takes the same bin values as the population, and iterates through
 * each sample index. The sample index is read from the raster band using the GDALRasterBand::RasterIO
 * function. The bins counts are updated accordingly. 
 *
 * @param RasterBandMetaData& band
 * @param std::vector<Index>* samples
 * @param std::vector<T> binVals
 * @param int nBins
 */
template <typename T>
std::vector<int64_t>
sampleDistribution(
	helper::RasterBandMetaData& band,
	std::vector<helper::Index>& samples,
	std::vector<T> binVals,
	int nBins)
{
	std::vector<int64_t> retval(nBins, 0);

	T val;
	for (const helper::Index& index : samples) {
		band.p_band->RasterIO(GF_Read, index.x, index.y, 1, 1, &val, 1, 1, band.type, 0, 0);
		for (int i = nBins - 1; i >= 0; i--) {
			if (binVals[i] <= val) {
				retval[i]++;
				break;
			}
		}
	}

	return retval;
}

/**
 * @ingroup dist
 *
 * This is the function which does most of calculations for the distribution calculation.
 * First the findMinMax() function is called to find the minimum and maximum pixel values
 * within the raster. Then the setBins() function is called to set the bin values
 * depending on the minumum and maximum values within the raster. Next, the populationDistribution()
 * function is ran to determine the population distribution within the calculated raster bins. Finally,
 * if the user provided a sample layer to compare against, the sampleDistribution() function is called
 * to determine the distribution within that particular sample.
 *
 * Both the findMinMax() and populationDistribution() function are called multiple times, each with
 * a different chunk (where the chunk count and chunk size is determined by the number of threads).
 * 
 * @param RasterBandMetaData& band
 * @param std::vector<Index>& sampled
 * @param int height
 * @param int width
 * @param int nBins
 * @param std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>& retval
 */
template <typename T>
void 
calculateDist(
	helper::RasterBandMetaData& band,
	std::vector<helper::Index>& sampled,
	int height,
	int width,
	int nBins,
	std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>& retval)
{
	//determine min and max values in the raster
	T min;
	T max;
	findMinMax<T>(band, width, height, min, max);

	//call the setBins function which sets the vector<double> bins to return to the user,
	//and the vector<T> bins to use to create the distribution. There is a seperate vector<T>
	//bins object so that while iterating through every pixel, they don't have to be cast to
	//type double to check.
	std::vector<double> dbins;
	std::vector<T> tbins;
	setBins<T>(min, max, nBins, band.type, dbins, tbins);

	//determine the population distribution of the raster band
	std::vector<int64_t> counts(nBins, 0);
	populationDistribution<T>(band, width, height, nBins, tbins, counts);	

	//add population distribution to return value
	retval.insert({std::string("population"), {dbins, counts}});

	//if samples are passed, add calculate sample distribution and add to return value
	if (sampled.size() != 0) {
		std::vector<int64_t> sampleCounts = sampleDistribution<T>(band, sampled, tbins, nBins);
		retval.insert({std::string("sample"), {dbins, sampleCounts}});
	}
}

/**
 * @ingroup dist
 * This function calculates the distribution of a particular raster band, and potentially
 * the distribution of a provided sample network of the same raster band.
 *
 * First, if there is a sample network provided, all of the points are iterated through
 * and converted to x/y index values usign the inverse geotransform (GDALInvGeoTransform()).
 *
 * After that, band metadata is acquired, and the calculateDistribution() function is called
 * with a template parameter depending on the raster type. After which, the return variable
 * (which was populated in calculateDistribution()) is returned. This variable is a map
 * where the key is either 'population' or 'sample' to specify whether the distribution of
 * the population or a sample distribution. The pair of vectors represent the bins and counts
 * per bin.
 *
 * @param GDALRasterWrapper *p_raster
 * @param int index
 * @param GDALVectorWrapper *p_vector
 * @param std::string layer
 * @param int nBins
 * @returns std::unordered_map<std::string, std::pari<std::vector<double>, std::vector<int64_t>>>
 */
std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>
dist(
	raster::GDALRasterWrapper *p_raster,
	int index,
	vector::GDALVectorWrapper *p_vector,
	std::string layer,
	int nBins)
{
	std::mutex datasetMutex;
	double *GT = p_raster->getGeotransform();
	double IGT[6];
       	GDALInvGeoTransform(GT, IGT);	
	
	//get samples as vector if indices, if samples are given
	std::vector<helper::Index> sampled;
	if (p_vector) {
		OGRLayer *p_layer = p_vector->getDataset()->GetLayerByName(layer.c_str());
	
		//check to ensure spatial reference system of raster and vector match	
		OGRSpatialReference rastSRS;
		rastSRS.importFromWkt(p_raster->getDataset()->GetProjectionRef());
		OGRSpatialReference *p_sampSRS = p_layer->GetSpatialRef();
		if (!rastSRS.IsSame(p_sampSRS)) {
			throw std::runtime_error("existing sample vector and raster do not have the same spatial reference system.");	
		}

		for (const auto& p_feature : *p_layer) {
			OGRGeometry *p_geometry = p_feature->GetGeometryRef();
			switch (wkbFlatten(p_geometry->getGeometryType())) {
				case OGRwkbGeometryType::wkbPoint: {
					OGRPoint *p_point = p_geometry->toPoint();
					double xCoord = p_point->getX();
					double yCoord = p_point->getY();
					sampled.push_back(helper::Index(
						static_cast<int>(IGT[0] + xCoord * IGT[1] + yCoord * IGT[2]),
						static_cast<int>(IGT[3] + xCoord * IGT[4] + yCoord * IGT[5])
					));
					break;
				}
				case OGRwkbGeometryType::wkbMultiPoint: {
					for (const auto& p_point : *p_geometry->toMultiPoint()) {
						double xCoord = p_point->getX();
						double yCoord = p_point->getY();
						sampled.push_back(helper::Index(
							static_cast<int>(IGT[0] + xCoord * IGT[1] + yCoord * IGT[2]),
							static_cast<int>(IGT[3] + xCoord * IGT[4] + yCoord * IGT[5])
						));
					}
					break;
				}
				default:
					throw std::runtime_error("encountered a geometry which was not a Point MultiPoint.");
			}
		}
	}

	std::mutex bandMutex;

	helper::RasterBandMetaData band;
	band.p_band = p_raster->getRasterBand(index);
	band.type = p_raster->getRasterBandType(index);
	band.size = p_raster->getRasterBandTypeSize(index);
	band.p_mutex = &datasetMutex;
	band.nan = band.p_band->GetNoDataValue();
	band.p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
	band.p_mutex = &bandMutex;

	int height = p_raster->getHeight();
	int width = p_raster->getWidth();

	std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>> retval;
	switch (band.type) {
		case GDT_Int8: 
			calculateDist<int8_t>(band, sampled, height, width, nBins, retval);
			break;
		case GDT_UInt16: 
			calculateDist<uint16_t>(band, sampled, height, width, nBins, retval);
			break;
		case GDT_Int16: 
			calculateDist<int16_t>(band, sampled, height, width, nBins, retval);
			break;
		case GDT_UInt32:
			calculateDist<uint32_t>(band, sampled, height, width, nBins, retval);
			break;
		case GDT_Int32:
			calculateDist<int32_t>(band, sampled, height, width, nBins, retval);
			break;
		case GDT_Float32:
			calculateDist<float>(band, sampled, height, width, nBins, retval);
			break;
		case GDT_Float64:
			calculateDist<double>(band, sampled, height, width, nBins, retval);	
			break;
		default:
			throw std::runtime_error("raster pixel data type not supported.");
	}

	return retval;
}

} //dist
} //sgs
