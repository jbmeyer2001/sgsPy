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

typedef clhs {

typedef oneapi::dal::homogen_table				DALHomogenTable;

struct Point {
	void *p_features;
	int x;
	int y;
}

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
	uint64_t usize;

	std::vector<std::vector<T>> quantiles;
	std::vector<std::vector<T>> corr;

	int64_t nFeat;

	xso::xoshiro_4x64_plus *p_rng = nullptr;
	uint64_t mask = 0;

	public:
	CLHSDataManager(int nFeat, xso::xoshiro_4x64_plus *p_rng) {
		this->nFeat = nFeat;
		this->points = 0;
		this->fi = 0;
		this->size = MILLION;
		this->points.resize(MILLION * nFeat);
		this->x.resize(MILLION);
		this->y.resize(MILLION);

		this->p_rng = p_rng;
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
	finalize(std::vector<std::vector<T>> corr, std::vector<std::vector<T>> quantiles) {
		this->corr = corr;
		this->quantiles = quantiles

		this->x.resize(this->size);
		this->y.resize(this->size);
		this->points.resize(this->size * nFeat);
		this->usize = static_cast<int64_t>(size);
	
		//use bit twiddling to fill the mask
		this->mask = static_cast<uint64_t>(size);
		this->mask--;
		this->mask |= this->mask >> 1;
		this->mask |= this->mask >> 2;
		this->mask |= this->mask >> 4;
		this->mask |= this->mask >> 8;
		this->mask |= this->mask >> 16;
		this->mask |= this->mask >> 32;
	}

	inline void
	randomPoint(Point& point) {
		int64_t index = ((*p_rng)() >> 11) & mask;

		while (index > usize) {
			index = ((*p_rng)() >> 11) & mask;
		}
	
		point.p_features = &points[index * nFeat];
		point.x = x[index];
		point.y = y[index];
	}
};

template <typename T>
inline void
readRaster(
	std::vector<RasterBandMetaData>& bands,
	CLHSDataManager<T>& clhs,
	Access& access,
	RandValController& rand,
	GDALDataType type,
	size_t size,
	int width,
	int height,
	int count,
	int nSamp)
{
	std::vector<std::vector<T>> probabilities(count);
	std::vector<std::vector<T>> quantiles(count);

	for (int i = 0; i < count; i++) {
		probabilities[i].resize(nSamp);
		quantiles[i].resize(nSamp);

		for (int j = 0; j < numSamples; j++) {
			probabilities[i][j] = static_cast<T>(j + 1) / static_cast<T>(nSamp);
		}
	}

	int xBlockSize = bands[0].xBlockSize;
	int yBlockSize = bands[0].yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (width + yBlockSize - 1) / yBlockSize;

	double eps = .001;
	
	std::vector<T> corrBuffer(count * xBlockSize * yBlockSize);
	std::vector<std::vector<T>> quantileBuffers(count);
	for (int i = 0; i < count; i++) {
		quantileBuffers[i].resize(xBlockSize * yBlockSize);
	}

	//create descriptor for correlation matrix streaming calculation with oneDAL
	const auto cor_desc = oneapi::dal::covariance::descriptor{}.set_result_options(dal::covariance::result_options::cor_matrix);
	oneapi::dal::covariance::partial_compute_result<> partial_result;

	//create tasks for quantile streaming calculation with MKL
	std::vector<VSLSSTaskPtr> quantileTasks;
	int status;
	MKL_INT quant_order_n = count;
	MKL_INT p = 1;
	MKL_INT n = xBlockSize * yBlockSize;
	MKL_INT nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;

	//MKL functions have different versions for single/double precision floating point data
	if (type == GDT_Float64) {
		for (int i = 0; i < count; i++) {
			//reinterpret cast the pointer for compiler reasons
			vsldSSNewTask(
				&quantileTasks[i], 
				&p, 
				&n, 
				&xstorage, 
				reinterpret_cast<double *>(probabilities[i].data()), 
				0, 
				0
			);
		}
	}
	else {
		for (int i = 0; i < count; i++) {
			//reinterpret cast the pointer for compiler reasons
			vslsSSNewTask(
				&quantileTasks[i],
				&p,
				&n,
				&xstorage,
				reinterpret_cast<float *>(probabilities[i].data()),
				0,
				0
			);
		}
	}

	bool calledEditStreamQuantiles = false;
	void *p_data = reinterpret_cast<void *>(corrBuffer.data());
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			//get block size
			int Valid, yValid;
			bands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);

			//read bands into memory
			for (int i = 0; i < count ; i++) {
				CPLErr err = bands[i].p_band->RasterIO(
					GF_Read,
					xBlock * xBlockSize,
					yBlock * yBlockSize,
					xValid,
					yValid,
					(void *)((size_t)p_data + i * size),
					xValid,
					yValid,
					size,
					size * static_cast<size_t>(count),
					size * static_cast<size_t>(count) * static_cast<size_t>(xBlockSize) 
				);

				if (err) {
					throw std::runtime_error("Error reading data from raster band.");
				}
			}

			//calculate rand vals
			rand.calculateRandValues();

			//iterate through pixels
			nFeatures = 0;
			for (int y = 0; y < yValid; y++) {
				int index = y * xBlockSize;
				for (int x = 0; x < xValid) {
					bool isNan = false;
					for (int b = 0; b < count; b++) {
						T val = corrBuffer[index * count + b];
						isNan = std::isnan(val) || val == bands[b].nan;

						if (isNan) {
							break;
						}	

						quantileBuffers[b][nFeatures] = val;
						corrBuffer[nFeatures * count + b] = val;
					}	

					if (!isNan) {
						nFeatures++;

						if (rand.next()) {
							clhs.addPoint(
								&corrBuffer[nFeatures * count],
								xBlock * xBlockSize + x,
								yBlock * yBlockSize + y		
							);
						}
					}
					index++;
				}
			}

			if (nFeatures == 0) {
				continue;
			}
			n = nFeatures; //tasks look to n for number of features

			//MKL functions have different versions for single/double precision floating point data
			if (type == GDT_Float64) {
				if (!calledEditStreamQuantiles) {
					for (int i = 0; i < count; i++) {
						//reinterpret cast the pointers for compiler reasons
						status = vsldSSEditStreamQuantiles(
							quantileTasks[i],
							&quant_order_n,
							reinterpret_cast<double *>(probabilities[i].data()),
							reinterpret_cast<double *>(quantiles[i].data()),
							&nparams,
							&eps
						);
					}
				}
				for (int i = 0; i < count; i++) {
					status = vsldSSCompute(
						quantileTasks[i],
						VSL_SS_STREAM_QUANTS,
						VSL_SS_METHOD_SQUANTS_ZW_FAST
					);
				}
			}
			else { //type == GDT_Float32
				if (!calledEditStreamQuantiles) {
					for (int i = 0; i < count; i++) {
						//reinterpret cast the pointers for compiler reasons
						status = vslsSSEditStreamQuantiles(
							quantileTasks[i],
							&quant_order_n,
							reinterpret_cast<float *>(probabilities[i].data()),
							reinterpret_cast<float *>(quantiles[i].data()),
							&nparams,
							&eps
						);
					}
				}
				for (int i = 0; i < count; i++) {
					status = vslsSSCompute(
						quantileTasks[i],
						VSL_SS_STREAM_QUANTS,
						VSL_SS_METHOD_SQUANTS_ZW_FAST
					);
				}

			}

			//update correlation matrix calculations
			DALHomogenTable table = DALHomogenTable(corrBuffer.data(), nFeatures, count, [](const double *){}, oneapi::dal::data_layout::row_major);
			partial_result = oneapi::dal::partial_train(cor_desc, partial_result, table); 
		}
	}

	//update clhs data manager with quantiles
	clhs.setQuantiles(quantiles);
	
	//calculate and update clhs data manager with correlation matrix
	auto result = oneapi::dal::finalize_compute(cor_desc, partial_result);
	auto correlation = result.get_cor_matrix();

	int64_t rows = correlation.get_row_count();
	int64_t cols = correlation.get_column_count();

	oneapi::dal::row_accessor<const T> acc {correlation};

	std::vector<std::vector<T>> correlation(count);
	for (int i = 0; i < count; i++) {
		correlation[i].resize(count);
		row = acc.pull{i, i + 1};

		for (int j = 0; j < count; j++) {
			correlation[i][j] = row[j];
		}
	}

	clhs.setCorrelation(correlation);
}

/**
 *
 */
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
clhs(
	GDALRasterWrapper *p_raster, 
	int numSamples,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string tempFolder,
	std::string filename)
{
	GDALAllRegister();

	int width = p_raster->getWidth();
	int height = p_raster->getHeight();
	int count = p_raster->getBandCount();
	double *GT = p_raster->getGeotransform();
	
	std::vector<double> xCoords, yCoords;

	std::vector<RasterBandMetaData> bands(p_raster->getBandCount());
	for (int i = 0; i < count; i++) {
		bands[i].p_band = p_raster->getRasterBand(i);
		bands[i].type = p_raster->getRasterBandType(i);
		bands[i].size = p_raster->getRasterBandTypeSize(i);
		bands[i].nan = bands[i].p_band->GetNoDataValue();
		bands[i].p_band->GetBlockSize(&bands[i].xBlockSize, &bands[i].yBlockSize);
	}

	//create output dataset before doing anything which will take a long time in case of failure.
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
	if (!p_driver) {
		throw std::runtime_error("unable to create output sample dataset driver.");
	}
	GDALDataset *p_samples = p_driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	if (!p_samples) {
		throw std::runtime_error("unable to create output dataset with driver.");
	}
	OGRLayer *p_layer = p_samples->CreateLayer("samples", nullptr, wkbPoint, nullptr);
	if (!p_layer) {
		throw std::runtime_error("unable to create output dataset layer.");
	}

	Access access(
		p_access,
		p_raster,
		layerName,
		buffInner,
		buffOuter,
		true,
		tempFolder,
		bands[0].xBlockSize,
		bands[0].yBlockSize
	);

	//fast random number generator using xoshiro256+
	//https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf
	xso::xoshiro_4x64_plus rng;
	multiplier = getProbabilityMultiplier(p_raster, 4, MILLION * 10, false, access.area);

	RandValController rand(bands[0].xBlockSize, bands[0].yBlockSize, multiplier, &rng);

	//get data type for all bands
	GDALDataType type = GDT_Float32;
	for (const RasterBandMetaData& band : bands) {
		if (band.type == GDT_Float64) {
			type = GDT_Float64;
			break;
		}
	}

	if (type == GDT_Float64) {	
		//create instance of data management class
		CLHSDataManager<double> clhs(numSamples, count);

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRaster<double>(bands, clhs, access, rand, type, sizeof(double), width, height, count, numSamples);
	}
	else { //type == GDT_Float32		
		//create instance of data management class
		CLHSDataManager<float> clhs(numSamples, count);

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRasterSP<float>(bands, clhs, access, rand, type, sizeof(float), width, height, count, numSamples);
	}

}

} // typedef clhs
