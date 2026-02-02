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
	T& max,
	std::mutex& mutex,
	std::condition_variable& cv,
	bool& finished) 
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

	mutex.lock();
	finished = true;
	mutex.unlock();
	cv.notify_all();
}

template <typename T>
std::vector<T>
setBins(
	T min,
	T max,
	int nBins,
	GDALDataType type,
	std::vector<double>& bins)
{
	//determine bucket values as doubles to display to the user.
	double step = ((double)max - (double)min) / ((double)nBins);
	double cur = (double)min;
	for (int i = 0; i <= nBins; i++) {
		bins.push_back(cur);
		cur += step;
	}

	//set bucket values as the data type used so less data type conversion is necessary in later code
	std::vector<T> retval(nBins);
	if (type == GDT_Float32 || type == GDT_Float64) {
		for (int i = 0; i < nBins; i++) {
			retval[i] = static_cast<T>(bins[i]);
		}
	}
	else {
		retval[0] = min;
		cur += step;
		for (int i = 1; i < nBins; i++) {
			retval[i] = std::ceil(cur);
			cur += step;
		}
	}

	return retval;
}

template <typename T>
void
populationDistribution(
	helper::RasterBandMetaData& band,
	int width,
	int yBlockStart,
	int yBlockEnd,
	int nBins,
	std::vector<T>& bucketVals,
	std::vector<int64_t>& bucketCounts,
	std::mutex& mutex,
	std::condition_variable& cv,
	bool& finished)
{
	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;

	T nan = static_cast<T>(band.nan);
	void *p_data = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);

	for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
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
							if (bucketVals[i] <= val) {
								bucketCounts[i]++;	
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

	mutex.lock();
	finished = true;
	mutex.unlock();
	cv.notify_all();
}

template <typename T>
std::vector<int64_t>
sampleDistribution(
	helper::RasterBandMetaData& band,
	std::vector<helper::Index>& samples,
	std::vector<T> bucketVals,
	int nBins)
{
	std::vector<int64_t> retval(nBins, 0);

	T val;
	for (const helper::Index& index : samples) {
		band.p_band->RasterIO(GF_Read, index.x, index.y, 1, 1, &val, 1, 1, band.type, 0, 0);
		for (int i = nBins - 1; i >= 0; i--) {
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
	int nBins,
	std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>& retval,
	int nThreads)
{
	pybind11::gil_scoped_acquire acquire;
	boost::asio::thread_pool pool(nThreads);

	//determine number of chunks and chunks size for each thread. A chunk is a group of blocks.
	int yBlocks = (height + band.yBlockSize - 1) / band.yBlockSize;
	int chunksize = std::max(1, yBlocks / nThreads);
	int nChunks = (yBlocks + chunksize - 1) / chunksize;

	std::mutex mutex;
	std::condition_variable cv;

	//don't use std::vector<bool> because can't reference individual bools from it
	bool *chunksFinished = reinterpret_cast<bool *>(VSIMalloc2(nChunks, sizeof(bool)));

	//determine min and max values in the raster
	std::vector<T> mins(nChunks);
	std::vector<T> maxs(nChunks);
	for (int chunk = 0; chunk < nChunks; chunk++) {
		int yStartBlock = chunk * chunksize;
		int yEndBlock = std::min(yStartBlock + chunksize, yBlocks);
		chunksFinished[chunk] = false;

		boost::asio::post(pool, std::bind(findMinMax<T>,
			std::ref(band),
			width,
			yStartBlock,
			yEndBlock,
			std::ref(mins[chunk]),
			std::ref(maxs[chunk]),
			std::ref(mutex),
			std::ref(cv),
			std::ref(chunksFinished[chunk])
		));
	}

	//using cv, wait for all of the threads to finish calculating their chunk before
	//calculating the final min and max values	
	std::unique_lock lock(mutex);
	bool allChunksFinished = true;
	for (int chunk = 0; chunk < nChunks; chunk++) {
		allChunksFinished &= chunksFinished[chunk];
	}

	while (!allChunksFinished) {
		cv.wait(lock);
		allChunksFinished = true;
		for (int chunk = 0; chunk < nChunks; chunk++) {
			allChunksFinished &= chunksFinished[chunk];
		}
	}

	//calculate the final min/max values
	T min = *std::min_element(mins.begin(), mins.end());
	T max = *std::max_element(maxs.begin(), maxs.end());

	//call the setBins function which sets the vector<double> bins to return to the user,
	//and the vector<T> bins to use to create the distribution. There is a seperate vector<T>
	//bins object so that while iterating through every pixel, they don't have to be cast to
	//type double to check.
	std::vector<double> dbins;
	std::vector<T> tbins = setBins<T>(min, max, nBins, band.type, dbins);

	//determine the population distribution of the raster band
	std::vector<std::vector<int64_t>> chunkCounts(nBins);
	for (int chunk = 0; chunk < nChunks; chunk++) {
		int yStartBlock = chunk * chunksize;
		int yEndBlock = std::min(yStartBlock + chunksize, yBlocks);

		std::vector<int64_t> counts(nBins, 0);
		chunkCounts[chunk] = counts;

		chunksFinished[chunk] = false;

		boost::asio::post(pool, std::bind(populationDistribution<T>,
			std::ref(band),
			width,
			yStartBlock,
			yEndBlock,
			nBins,
			std::ref(tbins),
			std::ref(chunkCounts[chunk]),
			std::ref(mutex),
			std::ref(cv),
			std::ref(chunksFinished[chunk])
		));	
	}

	//using cv, wait for all of the threads to finish calculating their chunk before
	//calculating the final counts in each bin
	allChunksFinished = true;
	for (int chunk = 0; chunk < nChunks; chunk++) {
		allChunksFinished &= chunksFinished[chunk];
	}

	while (!allChunksFinished) {
		cv.wait(lock);
		allChunksFinished = true;
		for (int chunk = 0; chunk < nChunks; chunk++) {
			allChunksFinished &= chunksFinished[chunk];
		}
	}

	//calculate the counts in each bin
	std::vector<int64_t> counts(nBins, 0);
	for (int i = 0; i < nBins; i++) {
		for (int j = 0; j < nChunks; j++) {
			counts[i] += chunkCounts[j][i];
		}
	}

	retval.insert({std::string("population"), {dbins, counts}});

	if (sampled.size() != 0) {
		std::vector<int64_t> sampleCounts = sampleDistribution<T>(band, sampled, tbins, nBins);
		retval.insert({std::string("sample"), {dbins, sampleCounts}});
	}

	pybind11::gil_scoped_release release;
}

std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>>
dist(
	raster::GDALRasterWrapper *p_raster,
	int index,
	vector::GDALVectorWrapper *p_vector,
	std::string layer,
	int nBins,
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
			calculateDist<int8_t>(band, sampled, height, width, nBins, retval, nThreads);
			break;
		case GDT_UInt16: 
			calculateDist<uint16_t>(band, sampled, height, width, nBins, retval, nThreads);
			break;
		case GDT_Int16: 
			calculateDist<int16_t>(band, sampled, height, width, nBins, retval, nThreads);
			break;
		case GDT_UInt32:
			calculateDist<uint32_t>(band, sampled, height, width, nBins, retval, nThreads);
			break;
		case GDT_Int32:
			calculateDist<int32_t>(band, sampled, height, width, nBins, retval, nThreads);
			break;
		case GDT_Float32:
			calculateDist<float>(band, sampled, height, width, nBins, retval, nThreads);
			break;
		case GDT_Float64:
			calculateDist<double>(band, sampled, height, width, nBins, retval, nThreads);	
			break;
		default:
			throw std::runtime_error("raster pixel data type not supported.");
	}

	return retval;
}

} //dist
} //sgs
