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

#include <condition_variable>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

namespace sgs {
namespace dist {

//no threading to start out
template <typename T>
void 
findMinMax(
	helper::RasterBandMetaData& band, 
	int width, 
	int yStartBlock,
	int yEndBlock,
	T& min,
	T& max) 
{
	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;

	min = std::numeric_limits<T>::max();
	max = std::numeric_limits<T>::min();
	T nan = static_cast<T>(band.nan);
	void *p_data = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);

	//calculate raster band minimum and maximum values
	for (int yBlock = yStartBlock; yBlock < yEndBlock; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			
			//read the band into memory
			rasterBandIO(band, p_data, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);

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

template <typename T>
std::vector<T>
setBuckets(
	T min,
	T max,
	int nBuckets,
	GDALDataType type,
	std::vector<double>& buckets)
{
	//determine bucket values as doubles to display to the user.
	double step = ((double)max - (double)min) / ((double)nBuckets);
	double cur = (double)min;
	for (int i = 0; i <= nBuckets; i++) {
		buckets.push_back(cur);
		cur += step;
	}

	//set bucket values as the data type used so less data type conversion is necessary in later code
	std::vector<T> retval(nBuckets);
	if (type == GDT_Float32 || type == GDT_Float64) {
		for (int i = 0; i < nBuckets; i++) {
			retval[i] = static_cast<T>(buckets[i]);
		}
	}
	else {
		retval[0] = min;
		cur += step;
		for (int i = 1; i < nBuckets; i++) {
			retval[i] = std::ceil(cur);
			cur += step;
		}
	}

	return retval;
}

template <typename T>
std::vector<int64_t>
populationDistribution(
	helper::RasterBandMetaData& band,
	int width,
	int yBlockStart,
	int yBlockEnd,
	int nBuckets,
	std::vector<T>& bucketVals,
	std::vector<int64_t>& bucketCounts)
{
	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;

	T nan = static_cast<T>(band.nan);
	void *p_data = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);

	for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			
			//read the band into memory
			rasterBandIO(band, p_data, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);

			for (int y = 0; y < yValid; y++) {
				int index = y * width;
				for (int x = 0; x < xValid; x++) {
					T val = reinterpret_cast<T *>(p_data)[index];
					
					if (val != nan && !std::isnan(val)) {
						for (int i = nBuckets - 1; i >= 0; i--) {
							if (bucketVals[i] <= val) {
								bucketCounts[i]++;	
								break;
							}
						}
					}	
				}
			}	
		}
	}

	VSIFree(p_data);
}

template <typename T>
std::vector<int64_t>
sampleDistribution(
	helper::RasterBandMetaData& band,
	std::vector<helper::Index>& samples,
	std::vector<T> bucketVals,
	int nBuckets)
{
	std::vector<int64_t> retval(nBuckets, 0);

	T val;
	for (const helper::Index& index : samples) {
		band.p_band->RasterIO(GF_Read, index.x, index.y, 1, 1, &val, 1, 1, band.type, 0, 0);
		for (int i = nBuckets - 1; i >= 0; i--) {
			if (bucketVals[i] <= val) {
				retval[i]++;
				break;
			}
		}
	}

	return retval;
}

template <typename T>
void 
calculateDist(
	helper::RasterBandMetaData& band,
	std::vector<helper::Index>& sampled,
	int height,
	int width,
	int nBuckets,
	std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>& retval,
	int nThreads)
{
	nThreads = 1; //for now testing
	
	//determine number of chunks and chunks size for each thread. A chunk is a group of blocks.
	int yBlocks = (height + band.yBlockSize - 1) / band.yBlockSize;
	int chunksize = yBlocks / nThreads;
	int nChunks = (yBlocks + chunksize - 1) / chunksize;
	
	//determine min and max values in the raster
	std::vector<T> mins;
	std::vector<T> maxs;
	for (int chunk = 0; chunk < nChunks; chunk++) {
		int yStartBlock = chunk * chunksize;
		int yEndBlock = std::min(yStartBlock + chunksize, yBlocks);

		T min, max;
		findMinMax<T>(band, width, yStartBlock, yEndBlock, min, max);

		//acquire mutex for adjusting these
		mins.push_back(min);
		maxs.push_back(max);
		//release mutex for adjusting these
	}
	T min = *std::min_element(mins.begin(), mins.end());
	T max = *std::max_element(maxs.begin(), maxs.end());

	//call the setBuckets function which sets the vector<double> buckets to return to the user,
	//and the vector<T> buckets to use to create the distribution. There is a seperate vector<T>
	//buckets object so that while iterating through every pixel, they don't have to be case to
	//type double to check.
	std::vector<double> dbuckets;
	std::vector<T> tbuckets = setBuckets<T>(min, max, nBuckets, band.type, dbuckets);

	//determine the population distribution of the raster band
	std::vector<int64_t> totalCounts;
	for (int chunk = 0; chunk < nChunks; chunk++) {
		int yStartBlock = chunk * chunksize;
		int yEndBlock = std::min(yStartBlock + chunksize, yBlocks);

		std::vector<int64_t> counts(nBuckets, 0);
		populationDistribution<T>(band, width, yStartBlock, yEndBlock, nBuckets, tbuckets, counts);

		//acquire mutex for adjusting these
		for (int i = 0; i < nBuckets; i++) {
			totalCounts[i] += counts[i];
		}
		//release mutex for adjusting these	
	}

	std::string key = std::string("population_") + std::string(band.p_band->GetDescription());
	retval.insert({key, {dbuckets, totalCounts}});

	if (sampled.size() != 0) {
		std::vector<int64_t> sampleCounts = sampleDistribution<T>(band, sampled, tbuckets, nBuckets);
		key = std::string("sample_") + std::string(band.p_band->GetDescription());
		retval.insert({key, {dbuckets, sampleCounts}});
	}
}

std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>
dist(
	raster::GDALRasterWrapper *p_raster,
	int index,
	vector::GDALVectorWrapper *p_vector,
	std::string layer,
	int nBuckets,
	int nThreads)
{
	std::mutex datasetMutex;
	double *GT = p_raster->getGeotransform();
	double IGT[6];
       	GDALInvGeoTransform(GT, IGT);	
	
	//get samples as vector if indices, if samples are given
	std::vector<helper::Index> sampled;
	if (p_vector) {
		OGRLayer *p_layer = p_vector->getDataset()->GetLayerByName(layer.c_str());
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

	helper::RasterBandMetaData band;
	band.p_band = p_raster->getRasterBand(index);
	band.type = p_raster->getRasterBandType(index);
	band.size = p_raster->getRasterBandTypeSize(index);
	band.p_mutex = &datasetMutex;
	band.nan = band.p_band->GetNoDataValue();
	band.p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
	
	int height = p_raster->getHeight();
	int width = p_raster->getWidth();

	std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>> retval;
	switch (band.type) {
		case GDT_Int8: 
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);
			break;
		case GDT_UInt16: 
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);
			break;
		case GDT_Int16: 
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);
			break;
		case GDT_UInt32:
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);
			break;
		case GDT_Int32:
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);
			break;
		case GDT_Float32:
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);
			break;
		case GDT_Float64:
			calculateDist<int8_t>(band, sampled, height, width, nBuckets, retval, nThreads);	
			break;
		default:
			throw std::runtime_error("raster pixel data type not supported.");
	}

	return retval;
}

} //dist
} //sgs
