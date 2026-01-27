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
 * @ingroup calculate
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
template typename<T>
findMinMax<T>(helper::RasterBandMetaData& band, int height, int width, int nBuckets) {
	int xBlocks = (width + band.xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	T tmin = std::numeric_limits<T>::max();
	T tmax = std::numeric_limits<T>::min();
	T nan = static_cast<T>(band.nan);
	p_data = VSIMalloc3(xBlockSize, yBlockSize, band.size);

	//calculate raster band minimum and maximum values
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			
			//read the band into memory
			rasterBandIO(band, p_data, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);

			for (int y = 0; y < yValid; y++) {
				int index = y * width;
				for (int x = 0; x < xValid; x++) {
					T val = p_data[index];
					if (val != nan && !std::isnan(val)) {
						tmin = std::min(tmin, val);
						tmax = std::max(tmax, val);
					}
					index++;
				}
			}	
		}
	}
	
	VSIFree(p_data);

	//assign bucket values using minimum, maximum, and total bucket count
	std::vector<T> retval(nBuckets + 1);
	if (band.type == GDT_Float32 || band.type == GDT_Float64) {
		T step = (tmax - tmin) / nBuckets;
		T cur = tmin;

		for (int i = 0; i < retval; i++) {
			retval[i] = cur;
			cur += step;
		}
	}
	else {
		min = static_cast<double>(tmin);
		max = static_cast<double>(tmax);
		double step = (max - min) / nBuckets;
		double cur = min;

		retval[0] = tmin;
		cur += step;
		for (int i = 1; i < retval; i++) {
			retval[i] = std::ceil(cur);
			cur += step;
		}
	}
	retval[nbuckets] = tmax;

	return retval;
}

template <typename T>
std::vector<int64_t>
populationDistribution(
	helper::RasterBandMetaData& band,
	int height,
	int width,
	std::vector<T> bucketVals,
	int nBuckets)
{
	int xBlocks = (width + band.xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	T nan = static_cast<T>(band.nan);
	T *p_data = VSIMalloc3(xBlockSize, YBlockSize, band.size);

	std::vector<int64_t> retval(nBuckets, 0);

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			
			//read the band into memory
			rasterBandIO(band, p_data, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);

			for (int y = 0; y < yValid; y++) {
				int index = y * width;
				for (int x = 0; x < xValid; x++) {
					T val = p_data[index];
					
					if (val != nan && !std::isnan(val)) {
						for (int i = nBuckets - 1; i >= 0; i--) {
							if (bucketVals[i] <= val) {
								retval[i]++;	
								break;
							}
						}
					}	
				}
			}	
		}
	}

	VSIFree(p_data);
	return retval;
}

template <typename T>
std::vector<int64_t>
sampleDistribution(
	helper::RasterBandMetaData& band,
	std::vector<Index>& samples,
	std::vector<T> bucketVals,
	int nBuckets)
{
	std::vector<int64_t> retval(nBuckets, 0);

	T val;
	for (const Index& index : samples) {
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


std::unordered_map<std::string, std::pair<std::vector<int>, std::vector<int>>>
dist(
	raster::GDALRasterWrapper *p_raster,
	std::vector<int> bandIndexes,
	vector::GDALVectorWrapper *p_vector,
	std::string layer,
	int nBuckets)
{
	std::mutex datasetMutex;
	double *GT = p_raster->getGeotransform();
	double IGT[6];
       	GDALInvGeoTransform(GT, IGT);	
	
	//get samples as vector if indices, if samples are given
	std::vector<Index> sampledIndexes;
	if (p_vector) {
		OGRLayer *p_layer = p_vector->getDataset()->GetLayerByName(layer.c_str());
		for (const auto& p_feature : *p_layer) {
			OGRGeometry *p_geometry = p_feature->GetGeometryRef();
			switch (wkbFlatter(p_geometry->getGeometryType())) {
				case OGRwkbGeometryType::wkbPoint: {
					OGRPoint *p_point = p_geometry->toPoint();
					double xCoord = p_point->GetX();
					double yCoord = p_point->GetY();
					sampledIndexes.push_back(Index(
						IGT[0] + xCoord * IGT[1] + yCoord * IGT[2],
						IGT[3] + xCoord * IGT[4] + yCoord * IGT[5]
					));
					break;
				}
				case OGRwkbGeometryType::skbMultiPoint: {
					for (const auto& p_point : *p_geometry->toMultiPoint()) {
						double xCoord = p_point->GetX();
						double yCoord = p_point->GetY();
						sampledIndexes.push_back(Index(
							IGT[0] + xCoord * IGT[1] + yCoord * IGT[2],
							IGT[3] + xCoord * IGT[4] + yCoord * IGT[5]
						));
					}
					break;
				}
				case default:
					throw std::runtime_error("encountered a geometry which was not a Point MultiPoint.");
			}
		}
	}

	std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<int64_t>>> retval;
	for (size_t i = 0; i < bandIndexes.size(); i++) {
		int index = bandIndexes[i];
		band.p_band = p_raster->GetRasterBand(index);
		band.type = p_raster->getRasterBandType(index);
		band.size = p_raster->getRasterBandTypeSize(index);
		band.p_mutex = &datasetMutex;
		band.nan = band.p_band->GetNoDataValue();
		band.p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
		
		switch (band.type) {
			case GDT_Int8: {
				std::vector<int8_t> buckets = findMinMax(band, height, width, nBuckets);
				std::vector<int64_t> counts = populationDistribution<int8_t>(band, height, width, buckets, nBuckets);
				
				break;
			}
			case GDT_UInt16: {
				
				break;
			}
			case GDT_Int16: {
				
				break;
			}
			case GDT_UInt32: {
				
				break;
			}
			case GDT_Int32: {
				
				break;
			}
			case GDT_Float32: {
				
				break;
			}
			case GDT_Float64: {
				
				break;
			}

		}
	}

}

} //dist
} //sgs
