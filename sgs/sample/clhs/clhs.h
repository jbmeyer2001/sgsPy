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

#include "utils/access.h"
#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

#include <mkl.h>
#include "oneapi/dal.hpp"
#include <xoshiro.h>

#define MILLION 1000000

typedef oneapi::dal::homogen_table				DALHomogenTable;

namespace clhs {

template <typename T>
struct Point {
	T *p_features = nullptr;
	int x = -1;
	int y = -1;
};

template <typename T>
inline size_t 
getQuantile(T val, std::vector<T>& quantiles) {
	auto it = std::lower_bound(quantiles.begin(), quantiles.end(), val);
	return (it == quantiles.end()) ?
		quantiles.size() :
		std::distance(quantiles.begin(), it);
}

template <typename T>
class CLHSDataManager {
	private:
	std::vector<T> features;
	std::vector<int> x;
	std::vector<int> y;
	size_t fi; //features index
	size_t count;
	int64_t size;
	uint64_t ucount;

	std::vector<std::vector<T>> quantiles;
	std::vector<std::vector<T>> corr;

	int nFeat;
	int nSamp;

	xso::xoshiro_4x64_plus *p_rng = nullptr;
	uint64_t mask = 0;

	public:
	CLHSDataManager(int nFeat, int nSamp, xso::xoshiro_4x64_plus *p_rng) {
		this->nFeat = nFeat;
		this->nSamp = nSamp;
		this->count = 0;
		this->fi = 0;
		this->size = MILLION;
		this->features.resize(MILLION * nFeat);
		this->x.resize(MILLION);
		this->y.resize(MILLION);

		this->p_rng = p_rng;
	}

	inline void
	addPoint(T *p_features, int x, int y) {
		for (int f = 0; f < nFeat; f++) {
			features[fi] = p_features[f];
			fi++;
		}

		this->x[this->count] = x;
	       	this->y[this->count] = y;
		this->count++;

		if (this->count == this->size) {
			this->features.resize(this->features.size() + MILLION * this->nFeat);
			this->x.resize(this->x.size() + MILLION);
			this->y.resize(this->y.size() + MILLION);
			this->size += MILLION;
		}	
	}

	inline void
	finalize(std::vector<std::vector<T>> corr) {
		if (this->count < static_cast<size_t>(this->nSamp)) {
			throw std::runtime_error("not enough points saved during raster iteration to conduct clhs sampling.");
		}
		
		this->corr = corr;

		this->x.resize(this->size);
		this->y.resize(this->size);
		this->features.resize(this->count * nFeat);
		this->ucount = static_cast<uint64_t>(this->count);

		//use bit twiddling to fill the mask
		this->mask = static_cast<uint64_t>(this->count);
		this->mask--;
		this->mask |= this->mask >> 1;
		this->mask |= this->mask >> 2;
		this->mask |= this->mask >> 4;
		this->mask |= this->mask >> 8;
		this->mask |= this->mask >> 16;
		this->mask |= this->mask >> 32;
	}

	inline uint64_t
	randomIndex() {
		uint64_t index = ((*p_rng)() >> 11) & mask;

		while (index > ucount) {
			index = ((*p_rng)() >> 11) & mask;
		}

		return index;
	}

	inline void
	getPoint(Point<T>& point, uint64_t index) {		
		point.p_features = this->features.data() + (index * nFeat);
		point.x = x[index];
		point.y = y[index];
	}

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
	std::vector<std::vector<T>>& quantiles,
	size_t size,
	int width,
	int height,
	int count,
	int nSamp)
{
	std::vector<std::vector<T>> probabilities(count);
	quantiles.resize(count);

	for (int i = 0; i < count; i++) {
		probabilities[i].resize(nSamp - 1);
		quantiles[i].resize(nSamp - 1);

		for (int j = 0; j < nSamp - 1; j++) {
			probabilities[i][j] = static_cast<T>(j + 1) / static_cast<T>(nSamp);
		}
	}

	int xBlockSize = bands[0].xBlockSize;
	int yBlockSize = bands[0].yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	double deps = .001;
	float seps = .001;

	std::vector<T> corrBuffer(count * xBlockSize * yBlockSize);
	std::vector<std::vector<T>> quantileBuffers(count);
	for (int i = 0; i < count; i++) {
		quantileBuffers[i].resize(xBlockSize * yBlockSize);
	}

	if (access.used) {
		access.band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, access.band.size); 
	}

	//create descriptor for correlation matrix streaming calculation with oneDAL
	const auto cor_desc = oneapi::dal::covariance::descriptor{}.set_result_options(oneapi::dal::covariance::result_options::cor_matrix);
	oneapi::dal::covariance::partial_compute_result<> partial_result;

	//create tasks for quantile streaming calculation with MKL
	std::vector<VSLSSTaskPtr> quantileTasks(count);
	int status;
	MKL_INT quant_order_n = quantiles[0].size();
	MKL_INT p = 1;
	MKL_INT n = xBlockSize * yBlockSize;
	MKL_INT nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;

	//MKL functions have different versions for single/double precision floating point data
	if (type == GDT_Float64) {
		for (int i = 0; i < count; i++) {
			//reinterpret cast the pointer for compiler reasons
			status = vsldSSNewTask(
				&quantileTasks[i], 
				&p, 
				&n, 
				&xstorage, 
				reinterpret_cast<double *>(quantileBuffers[i].data()), 
				0, 
				0
			);
		}
	}
	else {
		for (int i = 0; i < count; i++) {
			//reinterpret cast the pointer for compiler reasons
			status = vslsSSNewTask(
				&quantileTasks[i],
				&p,
				&n,
				&xstorage,
				reinterpret_cast<float *>(quantileBuffers[i].data()),
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
			int xValid, yValid;
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
					type,
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
				rasterBandIO(access.band, access.band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
			}

			//iterate through pixels
			n = 0;
			for (int y = 0; y < yValid; y++) {
				int index = y * xBlockSize;
				for (int x = 0; x < xValid; x++) {
					bool isNan = false;
					for (int b = 0; b < count; b++) {
						T val = corrBuffer[index * count + b];
						isNan = std::isnan(val) || val == bands[b].nan;

						if (isNan) {
							break;
						}	

						quantileBuffers[b][n] = val;
						corrBuffer[n * count + b] = val;
					}	

					if (!isNan) {
						n++;

						if ((!access.used || p_access[index] == 1) && rand.next()) {
							clhs.addPoint(
								corrBuffer.data() + (n * count),
								xBlock * xBlockSize + x,
								yBlock * yBlockSize + y		
							);
						}
					}
					index++;
				}
			}

			if (n == 0) {
				continue;
			}

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
							&deps
						);
					}
					calledEditStreamQuantiles = true;
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
							&seps
						);
					}
					calledEditStreamQuantiles = true;
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
			DALHomogenTable table = DALHomogenTable(corrBuffer.data(), n, count, [](const T *){}, oneapi::dal::data_layout::row_major);
			partial_result = oneapi::dal::partial_compute(cor_desc, partial_result, table); 
		}
	}

	if (access.used) {
		VSIFree(access.band.p_buffer);
	}

	//calculate and update clhs data manager with correlation matrix
	auto result = oneapi::dal::finalize_compute(cor_desc, partial_result);
	auto correlation = result.get_cor_matrix();

	oneapi::dal::row_accessor<const T> acc {correlation};

	n = 0;
	std::vector<std::vector<T>> corr(count);
	for (int i = 0; i < count; i++) {
		corr[i].resize(count);
		auto row = acc.pull({i, i + 1});

		for (int j = 0; j < count; j++) {
			corr[i][j] = row[j];
		}

		status = (type == GDT_Float64) ?
			vsldSSCompute(quantileTasks[i], VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZW) :
			vslsSSCompute(quantileTasks[i], VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZW);
		status = vslSSDeleteTask(&quantileTasks[i]);
	}

	clhs.finalize(corr);
}

template <typename T>
inline void
selectSamples(std::vector<std::vector<T>>& quantiles,
	      CLHSDataManager<T>& clhs,
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
	std::uniform_int_distribution<int> indexDist(0, nSamp);

	std::vector<std::vector<T>> corr(nFeat);
	std::vector<std::vector<int>> sampleCountPerQuantile(nFeat);
	std::vector<std::vector<std::unordered_set<uint64_t>>> samplesPerQuantile(nFeat);
	for (int i = 0; i < nFeat; i++) {
		sampleCountPerQuantile[i].resize(nSamp, 0);
		samplesPerQuantile[i].resize(nSamp);
		corr[i].resize(nFeat);
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
		
		if (!indicesMap.contains(index)) {
			indicesMap.emplace(index, i);
			
			Point<T> p;
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
	const auto cor_desc = oneapi::dal::covariance::descriptor{}.set_result_options(oneapi::dal::covariance::result_options::cor_matrix);		
	const auto result = oneapi::dal::compute(cor_desc, table);
	oneapi::dal::row_accessor<const T> acc {result.get_cor_matrix()};
	for (int i = 0; i < nFeat; i++) {
		auto row = acc.pull({i, i + 1});

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

	int test1 = 0;
	int test2 = 0;

	//begin annealing schedule. If we have a perfect latin hypercube -- or if we pass enough iterations -- stop iterating.
	while (temp > 0 && objQ != 0) {
		std::cout << "START" << std::endl;
		if (temp == 0.25 || temp == 0.50 || temp == 0.75 || temp == 1.0) {
			std::cout << "obj: " << obj << std::endl;
		}
		uint64_t swpIndex; //the index of the sample
		int i; //the index within the indices, x, y, and features vector so we know what to swap without searching
		if (dist(rng) < 0.5) {
			test1++;
			//50% of the time, choose a random sample to replace
			i = indexDist(rng);
			swpIndex = indices[i];
		}
		else {
			test2++;
			//50% of the time, choose the worst sample to replace
			int f, q, max;

			//get feature and quantile to remove
			for (f = 0; f < nFeat; f++) {	
				max = 0; 
				q = 0;

				for (int s = 0; s < nSamp; s++) {
					int count = sampleCountPerQuantile[f][s];

					if (count > max) {
						max = count;
						q = s;
					}
				}

				if (max != 1) {
					break;
				}
			}

			swpIndex = *samplesPerQuantile[f][q].begin();
			i = indicesMap.find(swpIndex)->second;
		}

		std::vector<T> oldf(nFeat);
		std::memcpy(
			reinterpret_cast<void *>(oldf.data()), 				//dst
			reinterpret_cast<void *>(features.data() + (i * nFeat)), 	//src
			nFeat * sizeof(T)						//size bytes
		);

		Point<T> p;
		uint64_t newIndex = clhs.randomIndex();
		while (indicesMap.contains(newIndex)) {
			newIndex = clhs.randomIndex();
		}
		clhs.getPoint(p, newIndex);
		std::memcpy(
			reinterpret_cast<void *>(features.data() + (i * nFeat)), 	//dst
			reinterpret_cast<void *>(p.p_features), 			//src
			nFeat * sizeof(T)						//size bytes
		);
			
		//recalculate sample count per quantile
		std::vector<int> oldq(nFeat);
		std::vector<int> newq(nFeat);
		for (int f = 0; f < nFeat; f++) {
			int q = getQuantile(oldf[f], quantiles[f]);
			oldq[f] = q;
			sampleCountPerQuantile[f][q]--;

			q = getQuantile(p.p_features[f], quantiles[f]);
			if (q >= nSamp) {
				std::cout << "q is larger than should be possible..." << std::endl;
				std::cout << "feature num is: " << f << std::endl;
				std::cout << "feature is: " << p.p_features[f] << std::endl;
				std::cout << "q: " << q << std::endl;
				std::cout << "quantiles: " << std::endl;
				for (int k = 0; k < nSamp; k++) {
					std::cout << "quantiles[" << k << "] = " << quantiles[f][k] << std::endl;
				}
			}
			newq[f] = q;
			sampleCountPerQuantile[f][q]++;
		}

		//recalculate objective function from quantiles
		T newObjQ = clhs.quantileObjectiveFunc(sampleCountPerQuantile);
		
		//recalculate correlation matrix
		const auto result = oneapi::dal::compute(cor_desc, table); // we update the table in place
		oneapi::dal::row_accessor<const T> acc {result.get_cor_matrix()};
		for (int j = 0; j < nFeat; j++) {
			auto row = acc.pull({j, j + 1});

			for (int k = 0; k < nFeat; k++) {
				corr[j][k] = row[k];
			}
		}

		//recalculate objective function from correlation matrix
		T newObjC = clhs.correlationObjectiveFunc(corr);

		T newObj = newObjQ + newObjC;
		T delta = newObj - obj;

		bool keep = dist(rng) < std::exp(-1 * delta * temp);

		if (keep) {
			std::cout << "KEEP" << std::endl;
			//update the new changes
			x[i] = p.x;
			y[i] = p.y;
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
			std::cout << "DONT KEEP" << std::endl;
			//revert anything already changed
			for (int f = 0; f < nFeat; f++) {
				sampleCountPerQuantile[f][newq[f]]--;
				sampleCountPerQuantile[f][oldq[f]]++;
			}

			std::memcpy(
				reinterpret_cast<void *>(features.data() + (i * nFeat)),
				reinterpret_cast<void *>(oldf.data()),
				nFeat * sizeof(T)
			);
		}

		//update annealing temperature
		temp -= d;
		std::cout << "TEMP: " << temp << std::endl;
	}

	std::cout << "OUT" << std::endl;
	//add samples to output layer
	for (int i = 0 ; i < nSamp; i++) {
		double xCoord = GT[0] + x[i] * GT[1] + y[i] * GT[2];
		double yCoord = GT[3] + x[i] * GT[4] + y[i] * GT[5];
		OGRPoint point = OGRPoint(xCoord, yCoord);
		addPoint(&point, p_layer);

		if (plot) {
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
		}
	}
}

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
	uint64_t multiplier = getProbabilityMultiplier(p_raster, 4, MILLION * 10, false, access.area);
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
		CLHSDataManager<double> clhs(nFeat, nSamp, &rng);

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRaster<double>(bands, clhs, access, rand, type, quantiles, sizeof(double), width, height, nFeat, nSamp);

		//select samples and add them to output layer
		selectSamples<double>(quantiles, clhs, rng, iterations, nSamp, nFeat, p_layer, GT, plot, xCoords, yCoords);
	}
	else { //type == GDT_Float32	
		std::vector<std::vector<float>> quantiles;

		//create instance of data management class
		CLHSDataManager<float> clhs(nFeat, nSamp, &rng);

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRaster<float>(bands, clhs, access, rand, type, quantiles, sizeof(float), width, height, nFeat, nSamp);

		//select samples and add them to output layer
		selectSamples<float>(quantiles, clhs, rng, iterations, nSamp, nFeat, p_layer, GT, plot, xCoords, yCoords);
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

	return {{xCoords, yCoords}, p_sampleVectorWrapper};
}

}
