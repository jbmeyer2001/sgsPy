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
	/**
	 *
	 */
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

	/**
	 *
	 */
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

	/**
	 *
	 */
	inline void
	finalize(std::vector<std::vector<T>> corr) {
		this->corr = corr;

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

	/**
	 *
	 */
	inline uint64_t
	randomIndex() {
		uint64_t index = ((*p_rng)() >> 11) & mask;

		while (index > usize) {
			index = ((*p_rng)() >> 11) & mask;
		}

		return index;
	}

	/**
	 *
	 */
	inline void
	getPoint(Point& point, uint64_t index) {		
		point.p_features = points.data() + (index * nFeat);
		point.x = x[index];
		point.y = y[index];
	}

	/**
	 *
	 */
	inline T
	quantileObjectiveFunc(std::vector<std::vector<int>>& sampleCountPerQuantile) {
		int retval = 0;

		for (const std::vector<int>& quantiles : sampleCountPerQuantile) {
			for (const int& count : quantiles) {
				retval += std::abs(count - 1);
			}
		}

		return static_cast<T>(retval);
	}

	/**
	 *
	 */
	inline T
	correlationObjectiveFunc(std::vector<std::vector<T>>& corr) {
		T retval = 0;

		for (size_t i = 0; i < this->corr.size(); i++) {
			for (size_t j = 0; j < this->corr[i].size(); j++) {
				retval += std::abs(corr[i][j] - this->corr[i][j]);
			}
		}
	
		return retval;
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
	std::vector<std::vector<T>>& quantiles;
	size_t size,
	int width,
	int height,
	int count,
	int nSamp)
{
	std::vector<std::vector<T>> probabilities(count);
	quantiles.resize(count);

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

	if (access.used) {
		access.band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, access.band.size); 
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

	int8_t *p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);

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

			//read access band into memory if used
			if (access.used) {
				rasterBandIO(access.band, access.band.p_buffer, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
			}

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

						if (access.used && p_access[index] == 1 && rand.next()) {
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
			DALHomogenTable table = DALHomogenTable(corrBuffer.data(), nFeatures, count, [](const T *){}, oneapi::dal::data_layout::row_major);
			partial_result = oneapi::dal::partial_train(cor_desc, partial_result, table); 
		}
	}

	if (access.used) {
		VSIFree(access.band.p_buffer);
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
		row = acc.pull({i, i + 1});

		for (int j = 0; j < count; j++) {
			correlation[i][j] = row[j];
		}
	}

	clhs.setCorrelation(correlation);
}

template <typename T>
inline void
selectSamples(std::vector<std::vector<T>>& quantiles,
	      CLHSDataManager& clhs,
	      xso::xoshiro_4x64_plus& rng,
	      int iterations,
	      int nSamp,
	      int nFeat,
	      OGRLayer *p_layer,
	      double *GT,
	      bool plot,
	      std::vector<double>& xCoords,
	      std::vector<double>& yCoords)
{
	std::uniform_real_distribution<T> dist(0.0, 1.0);
	std::unirofm_int_distribution<int> indexDist(0, nSamp);

	std::vector<std::vector<T>> corr(nFeat);
	std::vector<std::vector<T>> newCorr(nFeat);
	std::vector<std::vector<int>> sampleCountPerQuantile(nFeat);
	std::vector<std::vector<std::unordered_set<uint64_t>>> samplesPerQuantile(nFeat);
	for (int i = 0; i < nFeat; i++) {
		sampleCountPerQuantile[i].resize(nSamp, 0);
		samplesPerQuantile.resize(nSamp);
		corr.resize(nFeat);
		newCor.resize(nFeat);
	}

	std::vector<T> features(nSamp * nFeat);
	std::vector<int> x(nSamp);
	std::vector<int> y(nSamp);

	//indices vector used to get an index value in O(1) time after this vector is randomly indexed
	std::vector<uint64_t> indices(nSamp); 
	
	//indices map used to check whether an index is already used in O(1) time, and to keep track of it's index in the 'indices' vector
	std::unordered_map<uint64_t, int> indicesMap;

	//get first random samples
	int i = 0;
	while (i < nSamp) {
		uint64_t index = clhs.randomIndex();
		
		if (!indicesMap.find(index)) {
			indicesMap.emplace(index, i);
			
			Point p;
			clhs.getPoint(p, index);
		
			indices[i] = index;
			x[i] = p.x;
			y[i] = p.y;
			
			int fi = i * nFeat;
			for (int f = 0; f < nFeat; f++) {
				T val = p.p_features[f];
				features[fi + f] = val;

				int q = getQuantile<T>(val, quantiles[f]);
				sampleCountPerQuantile[f][q]++;
				samplesPerQuantile[f][q].insert(index);
			}

			i++;
		}
	}

	//define covariance calculation 
	DALHomogenTable table = DALHomogenTable(features.data(), nFeat, nSamp, [](const T *){}, oneapi::dal::data_layout::row_major);
	const auto cor_desc = oneapi::dal::covariance::descriptor{}.set_result_options(dal::covariance::result_options::cor_matrix);		
	const auto result = oneapi::dal::compute(cor_desc, table);
	oneapi::dal::row_accessor<const T> acc {result.get_cor_matrix()};
	for (int i = 0; i < nFeat; i++) {
		row = acc.pull({i, i + 1});

		for (int j = 0; j < nFeat; j++) {
			corr[i][j] = row[j];
		}
	}

	double temp = 1;
	double d = temp / static_cast<double>(iterations);

	T obj = 0;
	T objQ = clhs.quantileObjectiveFunc(sampleCountPerQuantile);
	T objC = clhs.correlationObjectiveFunc(corr);

	obj = objQ + objC;

	//begin annealing schedule. If we have a perfect latin hypercube -- or if we pass enough iterations -- stop iterating.
	while (t > 0 && objQ != 0) {
		uint64_t swpIndex; //the index of the sample
		int i; //the index within the indices, x, y, and features vector so we know what to swap without searching
		if (dist(rng) < 0.5) {
			//50% of the time, choose a random sample to replace
			i = indexDist(rng);
			swpIndex = indices[i];
		}
		else {
			//50% of the time, choose the worst sample to replace
			int feat = 0;

			//get feature and quantile to remove
			for (int f = 0; f < nFeat; f++) {	
				int max = 0; 
				int q = 0;

				for (int i = 0; i < nSamp; i++) {
					int count = sampleCountPerQuantile[f][i];

					if (count > max) {
						max = count;
						q = i;
					}
				}

				if (max != 1) {
					feat = f;
					break;
				}
			}

			swpIndex = *samplesPerQuantile[f][q].begin();
			i = indicesMap.find(swpIndex);
		}

		std::vector<T> oldf(nFeat);
		std::memcpy(
			reinterpret_cast<void *>(oldf.data()), 				//dst
			reinterpret_cast<void *>(features.data() + (i * nFeat)), 	//src
			nFeat * sizeof(T)						//size bytes
		);

		Point p;
		uint64_t newIndex = clhs.randomIndex();
		clhs.getPoint(newPoint, newIndex);
		std::memcpy(
			reinterpret_cast<void *>(features.data() + (i * nFeat)), 	//dst
			reinterpret_cast<void *>(p.p_features), 			//src
			nFeat * sizeof(T)						//size bytes
		);
			
		//recalculate sample count per quantile
		std::vector<int> oldq(nFeat);
		std::vector<int> newq(nFeat);
		for (int f = 0; f < nFeat; f++) {
			int q = getQuantile(oldf[f]);
			oldq[f] = q;
			sampleCountPerQuantile[f][q]--;

			q = getQuantile(p.p_features[f]);
			newq[f] = q;
			sampleCountPerQuantile[f][q]++;
		}

		//recalculate objective function from quantiles
		T newObjQ = clhs.quantileObjectiveFunc(sampleCountPerQuantile);
		
		//recalculate correlation matrix
		const auto result = oneapi::dal::compute(cor_desc, table); // we update the table in place
		oneapi::dal::row_accessor<const T> acc {result.get_cor_matrix()};
		for (int j = 0; j < nFeat; j++) {
			row = acc.pull({j, j + 1});

			for (int k = 0; k < nFeat; k++) {
				corr[j][k] = row[k];
			}
		}

		//recalculate objective function from correlation matrix
		T newObjC = clhs.correlationObjectiveFunc(newCorr);

		T newObj = newObjQ + newObjC;
		T delta = newObj - obj;

		bool keep = dist(rng) < std::exp(-1 * delta * temp);

		if (keep) {
			//update the new changes
			x[i] = point.x;
			y[i] = point.y;
			indicesMap.erase(indices[i]);
			indicesMap.emplace(newIndex, i);
				
			for (int f = 0; f < nFeat; f++) {
				samplesPerQuantile[f][oldq[f]].erase(indices[i]);
				samplesPerQuantile[f][newq[f]].insert(newIndex);
			}

			indices[i] = newIndex;
			
			objC = newObjC;
			objQ = newObjQ;
			obj = newObj;
		}
		else {
			//revert anything already changed
			for (int f = 0; f < nFeat; f++) {
				sampleCountPerQuantile[f][newq[f]]--;
				sampleCountPerquantile[f][oldq[f]]++;
			}

			std::memcpy(
				reinterpret_cast<void *>(features.data() + (i * nFeat)),
				reinterpret_cast<void *>(oldf.data()),
				nFeat * sizeof(T),
			);
		}

		//update annealing temperature
		temp -= d;
	}

	//add samples to output layer
	for (int i = 0 ; < nSamp; i++) {
		double x = GT[0] + x[i] * GT[1] + y[i] * GT[2];
		double y = GT[3] + x[i] * GT[4] + y[i] * GT[5];
		OGRPoint point = OGRPoint(x, y);
		addPoint(&point, p_layer);

		if (plot) {
			xCoords.push_back(x[i]);
			yCoords.push_back(y[i]);
		}
	}
}

/**
 *
 */
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *>
clhs(
	GDALRasterWrapper *p_raster, 
	int nSamp,
	int iterations,
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
	int nFeat = p_raster->getBandCount();
	double *GT = p_raster->getGeotransform();
	
	std::vector<double> xCoords, yCoords;

	std::vector<RasterBandMetaData> bands(p_raster->getBandCount());
	for (int i = 0; i < nFeat; i++) {
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
		std::vector<std::vector<double>> quantiles;
		
		//create instance of data management class
		CLHSDataManager<double> clhs(nSamp, nFeat);

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRaster<double>(bands, clhs, access, rand, quantiles, type, sizeof(double), width, height, nFeat, nSamp);

		//select samples and add them to output layer
		selectSamples<double>(quantiles, clhs, rng, iterations, nSamp, nFeat, p_layer, GT, plot, xCoords, yCoords);
	}
	else { //type == GDT_Float32		
		std::vector<std::vector<float>> quantiles;

		//create instance of data management class
		CLHSDataManager<float> clhs(nSamp, nFeat);

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRasterSP<float>(bands, clhs, access, rand, quantiles type, sizeof(float), width, height, nFeat, nSamp);

		//select samples and add them to output layer
		selectSamples<double>(quantiles, clhs, rng, iterations, nSamp, nFeat, p_layer, GT, plot, xCoords, yCoords);
	}

	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_samples);

	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	return {{sCoords, yCoords}, p_sampleVectorWrapper};
}

} // typedef clhs
