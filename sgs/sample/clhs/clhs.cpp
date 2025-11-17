/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of CLHS
 * Author: Joseph Meyer
 * Date: November, 2025
 *
 ******************************************************************************/

#include <iostream>
#include <random>

#include "access.h"
#include "helper.h"
#include "raster.h"
#include "vector.h"

#include <mkl.h>
#include "oneapi/dal.hpp"
#include <xoshiro.h>

#define MILLION 1000000

typedef oneapi::dal::homogen_table				DALHomogenTable;

/**
 *
 */
template <typename T>
inline size_t 
getQuantile(T val, std::vector<T>& quantiles) {
	auto it = std::lower_bound(quantiles.begin(), quantiles.end(), val);
	return (it == quantiles.end()) ?
		quantiles.size() :
		std::distance(quantiles.begin(), it);
}

/**
 *
 */
template <typename T>
class CLHSDataManager {
	private:
	std::vector<T> features;
	std::vector<int> x;
	std::vector<int> y;
	int64_t fi; //features index
	int64_t points;
	int64_t size;

	std::vector<T> samplePoints;
	std::vector<int> sampleX;
	std::vector<int> sampleY;

	std::vector<std::vector<int>> quantileCounts;
	std::vector<std::vector<T>> quantiles;
	std::vector<T> corr;

	int64_t nSamp;
	int64_t nFeat;

	public:
	CLHSDataManager(int nSamp, int nFeat) {
		this->nSamp = nSamp;
		this->nFeat = nFeat;
		this->points = 0;
		this->fi = 0;
		this->size = MILLION;

		this->quantileCounts.resize(nFeat);
		for (int f = 0; f < nFeat; f++) {
			this->quantileCounts[f].resize(nSamp, 0);
		}
		
		this->samplePoints.resize(nSamp);
		this->sampleX.resize(nSamp);
		this->sampleY.resize(nSamp);

		this->points.resize(MILLION * nFeat);
		this->x.resize(MILLION);
		this->y.resize(MILLION);	
	}

	inline void
	addPoint(T *p_features, int x, int y) {
		for (int64_t f = 0; f < nFeat; f++) {
			points[fi] = p_features[f];
			fi++;
		}

		x[points] = x;
	       	y[points] = y;
		points++;

		if (points == size) {
			points.resize(points.size + MILLION * nFeat);
			x.resize(x.size() + MILLION);
			y.resize(y.size() + MILLION);
			size += MILLION;
		}	
	}

	inline void
	setCorrelationMatrix(std::vector<T> corr) {
		this->corr = corr;
	}

	inline void
	setQuantiles(std::vector<T> quantiles) {
		this->quantiles = quantiles;
	}
};



/**
 *
 */
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
clhs(
	GDALRasterWrapper *p_raster, 
	int numSamples

)
