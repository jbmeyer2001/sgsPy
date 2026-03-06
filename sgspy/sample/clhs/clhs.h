/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of CLHS
 * Author: Joseph Meyer
 * Date: November, 2025
 *
 ******************************************************************************/

/**
 * @defgroup clhs clhs
 * @ingroup sample
 */

#include <iostream>
#include <random>

#include "utils/access.h"
#include "utils/existing.h"
#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

#include <mkl.h>
#include "oneapi/dal.hpp"
#include <xoshiro.h>

#define MILLION 1000000

typedef oneapi::dal::homogen_table				DALHomogenTable;

namespace sgs {
namespace clhs {

/**
 * @ingroup clhs
 * Structure for containing the x and y positions of a point,
 * along with an array of its feature values.
 */
template <typename T>
struct Point {
	T *p_features = nullptr;
	int x = -1;
	int y = -1;
};

/**
 * @ingroup clhs
 * This function is used to get the quantile which a particular value fits into.
 * It requires both that value, and a vector of the quantile values to be passed.
 *
 * @param T val
 * @param std::vector<T>& quantiles
 * @returns size_t
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
 * @ingroup clhs
 * This class is responsible for managing the data for the clhs sampling method.
 *
 * It contains a vector with the feature values of all added pixels, as well as 
 * the x and y values of those pixels, and can randomly return one of those pixels
 * as desired. It also stores the correlation matrix of the raster. 
 */
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

	//for existing sample points
	std::vector<T> efeatures;
	std::vector<int> ex;
	std::vector<int> ey;
	size_t efi = 0;
	size_t ecount = 0;

	std::vector<std::vector<T>> corr;

	int nFeat;
	int nSamp;

	xso::xoshiro_4x64_plus *p_rng = nullptr;
	uint64_t mask = 0;

	public:
	/**
	 * Constructor, sets the nFeat, nSamp, and random number generator. Also
	 * sizes the vectors to 1,000,000 points (initially). The vector will
	 * resize as required, and sizes down once raster reading is completed.
	 *
	 * set sizes and count values of existing pixels if required.
	 *
	 * @param int nFeat
	 * @param int nSamp
	 * @param xso::xoshiro_4x64_plus *p_rng
	 * @param int existingCount
	 */
	CLHSDataManager(int nFeat, int nSamp, xso::xoshiro_4x64_plus *p_rng, size_t existingCount) {
		this->nFeat = nFeat;
		this->nSamp = nSamp;
		this->count = 0;
		this->fi = 0;
		this->size = MILLION;
		this->features.resize(MILLION * nFeat);
		this->x.resize(MILLION);
		this->y.resize(MILLION);

		this->p_rng = p_rng;

		if (existingCount != 0) {
			this->efeatures.resize(existingCount * static_cast<size_t>(nFeat));
			this->ex.resize(existingCount);
			this->ey.resize(existingCount);
		}
	}

	/**
	 * This function saves a point to the data manager. It takes a feature
	 * vector (array), as a parameter alongside the x and y indices of the
	 * pixel containing those features.
	 *
	 * The x, y, and features vectors are updated accordingly and resized
	 * if required.
	 *
	 * @param T *p_features
	 * @param int x
	 * @param int y
	 */
	inline void
	addPoint(T *p_features, int x, int y) {
		for (int f = 0; f < nFeat; f++) {
			features[this->fi] = p_features[f];
			this->fi++;
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

	/**
	 * This function saves specifically an existing sample point to the data manager. It
	 * takes a feature vector (array) as a parameter alongside the x and y indices of the pixel
	 * containing those features.
	 *
	 * The x, y, and features vectors are updated accordingly.
	 *
	 * @param T* p_features
	 * @param int x
	 * @param int y
	 */
	inline void
	addExistingPoint(T *p_features, int x, int y) {
		for (int f = 0; f < nFeat; f++) {
			this->efeatures[this->efi] = p_features[f];
			this->efi++;
		}

		this->ex[this->ecount] = x;
		this->ey[this->ecount] = y;
		this->ecount++;
	}

	/**
	 * This function is called once the raster has been read, meaining no
	 * more points will be added to the data manager. The correlation
	 * matrix, which is calculated just after the raster reading finishes,
	 * is passed as a parameter so that it can be saved.
	 *
	 * The x, y, and features vectors are resized.
	 *
	 * A mask, which is used along with the random number generator to 
	 * generate random indices within the saved points, is generated.
	 * The mask value is all 1 starting from the most significant
	 * bit which is 1 when the index is at it's largest.
	 *
	 * When anded against a new random number, the mask value will generate
	 * a number which can be any of the indices, and is quite likely not to
	 * be larger than the capacity. If it is larger than the capacity it can
	 * just be calculated again.
	 *
	 * @param std::vector<std::vector<T>>& corr 
	 */
	inline void
	finalize(std::vector<std::vector<T>>& corr) {
		if (this->count < static_cast<size_t>(this->nSamp)) {
			throw std::runtime_error("not enough points saved during raster iteration to conduct clhs sampling.");
		}
		
		this->corr = corr;

		this->x.resize(this->count);
		this->y.resize(this->count);
		this->features.resize(this->count * nFeat);
		this->ucount = static_cast<uint64_t>(this->count);

		//use bit twiddling to fill the mask
		this->mask = static_cast<uint64_t>(this->count);
		this->mask |= this->mask >> 1;
		this->mask |= this->mask >> 2;
		this->mask |= this->mask >> 4;
		this->mask |= this->mask >> 8;
		this->mask |= this->mask >> 16;
		this->mask |= this->mask >> 32;

		//resize existing sample vectors 
		this->ex.resize(this->ecount);
		this->ey.resize(this->ecount);
		this->features.resize(this->ecount * nFeat);
	}

	/**
	 * This function is for generating a random index among the points saved.
	 *
	 * Due to using the random number generator and a mask direcly (faster
	 * than something like an std::uniform_int_distribution) there may
	 * be some values which are larger than the total number of indices
	 * which are occupied. Due to this, indices are generated until
	 * there's one which is a valid index.
	 *
	 * The reason the generator is initially bit shifted by 11 is because
	 * the type of generator used is not as random in the first 11 bits.
	 *
	 * @returns uint64_t
	 */
	inline uint64_t
	randomIndex() {
		uint64_t index = ((*p_rng)() >> 11) & mask;

		while (index >= this->ucount) {
			index = ((*p_rng)() >> 11) & mask;
		}

		return index;
	}

	/**
	 * This function sets the x, y, and features pointer for a given index.
	 *
	 * @param Point<T>& point
	 * @param uint64_t index
	 */
	inline void
	getRandomPoint(Point<T>& point) {
		uint64_t index = randomIndex();

		point.p_features = this->features.data() + (index * nFeat);
		point.x = x[index];
		point.y = y[index];
	}

	/**
	 * This function calculates the output to the continuous objective function
	 * for the clhs method. The goal of the clhs method is to have a latin hypercube
	 * output of samples. This is defined as having a single sample between each quantile
	 * of the hypercube for every feature. This function returns 0 if the sample is a 
	 * latin hypercube, and the number of samples off if it is not. 
	 *
	 * The output of this function, along with the output of the correlation objective
	 * function are used to determine whether a sample should be kept or not. Lower values
	 * are better.
	 *
	 * @param std::vector<std::vector<int>>& sampleCountPerQuantile
	 * @returns T
	 */
	inline T
	quantileObjectiveFunc(std::vector<int>& sampleCountPerQuantile) {
		int retval = 0;

		for (const int& count : sampleCountPerQuantile) {
			retval += std::abs(count - 1);
		}

		return static_cast<T>(retval);
	}

	/**
	 * This function calculates the output to the objective function for the 
	 * differences in the correlation matricies for the clhs method. This
	 * function penalizes the sample for having features which are correlated
	 * differently than the whole sample space.
	 *
	 * The output of this function, along with the output of the quantile
	 * objective function are used to determine whether as ample should be kept or not.
	 * Lower values are better.
	 *
	 * @param std::vector<std::vector<T>> corr
	 * @returns T
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

	/**
	 * Copy features, x, and y vectors from CLHS data manager.
	 *
	 * @param std::vector<T>& features
	 * @param std::vector<int>& x
	 * @param std::vector<int>& y
	 */
	inline void
	getExistingFeatures(std::vector<T>& features, std::vector<int>& x, std::vector<int>& y) {
		features.resize(this->efeatures.size());
		x.resize(this->ex.size());
		y.resize(this->ey.size());

		std::memcpy(features.data(), this->efeatures.data(), this->efeatures.size() * sizeof(T));
		std::memcpy(x.data(), this->ex.data(), this->ex.size() * sizeof(int));
		std::memcpy(y.data(), this->ey.data(), this->ey.size() * sizeof(int));
	}
};

/**
 * @ingroup clhs
 * This function is responsible for three different tasks as it reads through the raster.
 * The raster is read in blocks, so as to be as memory efficient as possible and 
 * avoid the potential issue of a raster which is too large to fit in memory.
 *
 * The three tasks are as follows:
 *  - calculate the quantile values for every feature, where the number of equally sized
 *  quantiles is the sample size.
 *  - calculate the correltion matrix of the features in the raster.
 *  - Save points, along with their x/y coordinates.
 *
 * QUANTILES:
 * MKL (Math Kernel Library) is used to calculate the quantiles using a quantile streaming
 * algorithm. As the pixels are iterated through, each pixel is a row an each feature
 * is a column -- this is because it is how the correlation matrix calculation expects
 * the values to be, but it also simplifies the removal of nan pixels. However,
 * the quantile calculation is done per-feature, so each feature has it's own vector of 
 * values which are copied over from the original matrix while reading. These are
 * in the quantileBuffers vector of vectors. 
 *
 * Additionally, there are seperate MKL functions for single precision vs double precision
 * floating point values. For this reason, the type is passed (in addition to being templated)
 * so that the correct MKL function is called according to the data type of the raster.
 *
 * Due to the template function also being passed, the pointer type passed to those type-specific
 * MKL functions must be casted to the correct type. This is because, the C++ compiler is not
 * aware that the float template will never be passed along with a type of GDT_Float64, and double
 * will never be passed along with a type of GDT_Float32.
 *
 * CORRELATION:
 * oneDAL is used to calculate the correlation matrix of the raster features using a streaming
 * algorithm. The features are read such that a pixel corresponds to a row and a feature
 * corresponds to a column in the buffer passed to the quantile streaming algorithm. This is to
 * simplify the removal of nan pixels, and to maximize cache-line efficiency.
 *
 * POINTS:
 * The most obvious way to calculate clhs would be to read the entire raster into memory, remove
 * the nan pixels, then run the simulated annealing algorithm over the entire raster. Essentially,
 * every pixel in the raster has an equal chance of being selected. This method however falls apart
 * when there are enough pixels so as to be encredibly memory-inefficent, or so many that they can't
 * fit in memory at once. In this case, we cannot store the features of every single pixel. 
 *
 * Instead, as we iterate through the raster we randomly select pixels which will be added to a pool. The pixels
 * within that pool are options for being added to the latin hypercube. The probability of each pixel
 * being added to this pool depends on the total size of the raster (as well as the number of accessible
 * pixels, if the access parameter is given), but every accessible pixel has an equal chance of being
 * included in the pool. The size goal for the pool is 10,000,000 pixels because it is (hopefully) enough
 * to give the simulated annealing algorithm all the options for quantiles it needs, without breaking memory.
 * If the raster has a small enough number of pixels however, they will all be added to the pool.
 *
 * @param std::vector<RasterBandMetaData>& bands
 * @param CLHSDataManager<T>& clhs
 * @param Access& access
 * @param RandValController& rand
 * @param GDALDataType type
 * @param std::vector<std::vector<T>>& quantiles
 * @param size_t size
 * @param int width
 * @param int height
 * @param int count
 * @param int nSamp
 */
template <typename T>
inline void
readRaster(
	std::vector<helper::RasterBandMetaData>& bands,
	CLHSDataManager<T>& clhs,
	access::Access& access,
	existing::Existing& existing,
	helper::RandValController& rand,
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
	float seps = .001f;

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
					T *p_buff = corrBuffer.data() + (index * count);
					for (int b = 0; b < count; b++) {
						T val = p_buff[b];
						isNan = std::isnan(val) || val == bands[b].nan;

						if (isNan) {
							break;
						}	

						quantileBuffers[b][n] = val;
						corrBuffer[n * count + b] = val;
					}	

					if (!isNan) {
						n++;
						
						bool accessible = !access.used || p_access[index] != 1;

						if (existing.used && existing.containsIndex(x + xBlock * xBlockSize, y + yBlock * yBlockSize)) {
							clhs.addExistingPoint(
								p_buff,
								xBlock * xBlockSize + x,
								yBlock * yBlockSize + y
							);
						}
						else if (accessible && rand.next()) {
							clhs.addPoint(
								p_buff,
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

/**
 * @ingroup clhs
 * This function is responsible for selecting the samples will be in the final sample
 * using the Conditioned Latin Hypercube Sampling (clhs) method. The CLHSDataManager
 * contains a pool of points which may be selected, as well as objective functions which
 * will calculate how good each sample is depending on the number of samples within
 * each features quantiles, and how closely matching the correlation matrices are between
 * the sample compared to the entire sample space.
 *
 * The clhs method uses a simulated annealing algorithm, wherein the sample is constanty being
 * updated and tested by replacing one of the samples. The old vs new samples are compared usign
 * the objective functions, and the likelyhood that the sample is added (whether it's an improvement
 * or not) depends both on a random probability, and the tempreature parameter which decreases after
 * every iteration. Essentially, at the start of the algorithm it is far more likely that a sample
 * may be accepted if it makes the objective function worse. The reason why the algorithm doesn't always
 * accept an improvement and reject a decrease is to avoid finding a local maximum (which is not a
 * global maximum). At the start of each iteration, there is also a 50% chance that the sample removed
 * is NOT random, but removed from the quantile which has the most samples (essentially removing the 
 * 'worst' pixel).
 *
 * Once all of the samples are chosen, they are added to the output layer and the function completes.
 *
 * @param std::vector<std::vector<T>>& quantiles
 * @param xso::xoshiro_4x64_plus& rng
 * @param existing::Existing& existing
 * @param int replace
 * @param int iterations
 * @param int nSamp
 * @param int nFeat
 * @param OGRLayer *p_layer
 * @param double *GT
 * @param bool plot
 * @param std::vector<double>& xCoords
 * @param std::Vector<double>& yCoords
 */
template <typename T>
inline void
selectSamples(std::vector<std::vector<T>>& quantiles,
	      CLHSDataManager<T>& clhs,
	      xso::xoshiro_4x64_plus& rng,
	      existing::Existing& existing,
	      size_t replace,
	      size_t iterations,
	      size_t nSamp,
	      size_t nFeat,
	      OGRLayer *p_layer,
	      double *GT,
	      bool plot,
	      std::vector<double>& xCoords,
	      std::vector<double>& yCoords)
{
	std::vector<T> features;
	std::vector<int> x;
	std::vector<int> y;

	std::vector<int> sampleCountPerQuantile(nFeat * nSamp, 0); //nFeat x nSamp 2d array
	std::vector<int> quantilesOfEachSample(nSamp * nFeat, 0); //nSamp x nFeat 2d array
	
	//if there are existing samples, add all of them.
	if (existing.used) {
		clhs.getExistingFeatures(features, x, y);
	}

	//Then, remove up to 'replace' number of them which are redundant
	//'redundant' samples are samples which cause over-representation in feature quantiles
	
	size_t neSamples = x.size(); //number of existing samples
	if (neSamples > 0 && replace != 0) {
		for (size_t si = 0; si < neSamples; si++) {
			for (size_t fi = 0; fi < nFeat; fi++) {
				T val = features[(si * nFeat) + fi];
				size_t q = getQuantile<T>(val, quantiles[fi]);
				sampleCountPerQuantile[(fi * nSamp) + q]++;
				quantilesOfEachSample[(si * nFeat) + fi] = q;
			}
		}

		while (replace > 0 && neSamples > 0) {
			size_t worstRedundancy = 0;
			size_t worstRedundancyIndex = 0;

			//get the sample with the worst redundancy. In other words, the one which is in the
			//most over-represented quantile across all features
			for (size_t si = 0; si < neSamples; si++) {
				size_t curSampleRedundancy = 0;
				for (size_t fi = 0; fi < nFeat; fi++) {
					size_t q = quantilesOfEachSample[(si * nFeat) + fi];
					curSampleRedundancy += sampleCountPerQuantile[(fi * nSamp) + q];
				}
				if (curSampleRedundancy > worstRedundancy) {
					worstRedundancy = curSampleRedundancy;
					worstRedundancyIndex = si;
				}	
			}

			//if the remaining samples are already forming a partial latin hypercube
			//(no more than 1 sample per quantile of each feature)
			//then don't remove any more indices!
			if (worstRedundancy == nFeat) {
				break;
			}

			//replace the most redundant sample with the last sample in the vector, and decrease
			//the size of the vector by 1. But first, adjust the values of sampleCountPerQuantile
			size_t si = worstRedundancyIndex;
			for (size_t fi = 0; fi < nFeat; fi++) {
				size_t q = quantilesOfEachSample[(si * nFeat) + fi];
				sampleCountPerQuantile[(fi * nSamp) + q]--;
			}

			x[si] = x[neSamples - 1];
			y[si] = y[neSamples - 1];
			std::memcpy(features.data() + (si * nFeat), 
				    features.data() + ((neSamples - 1) * nFeat), 
				    sizeof(T) * nFeat);
			std::memcpy(quantilesOfEachSample.data() + (si * nFeat),
				    quantilesOfEachSample.data() + ((neSamples - 1) * nFeat),
				    sizeof(int) * nFeat);

			neSamples--;
			replace--;
		}		
	}

	//NOW, the neSamples variable contains the number of existing samples which MUST be kept.
	//This number also happens to be the first index of the remaining space in the vector which
	//may be filled in with non-existing samples. If there were no existing samples this value
	//is 0.
	size_t starti = neSamples;
	std::unordered_set<uint64_t> points;

	//Add all of the existing samples to the output layer, and add to indices map
	helper::Field fieldExistingTrue("existing", 1);
	for (size_t si = 0; si < neSamples; si++) {
		OGRPoint point = existing.getPoint(x[si], y[si]);
		helper::addPoint(&point, p_layer, &fieldExistingTrue);

		if (plot) {
			xCoords.push_back(point.getX());
			yCoords.push_back(point.getY());
		}

		points.insert((((uint64_t) x[si]) << 32) | ((uint64_t) y[si]));
	}

	//if there are already enough (existing) samples, return and don't add any more
	if (starti >= nSamp) {
		return;
	}

	std::uniform_real_distribution<T> dist(0.0, 1.0);
	std::uniform_int_distribution<size_t> indexDist(starti, nSamp - 1);

	std::vector<std::vector<T>> corr(nFeat);
	for (int i = 0; i < nFeat; i++) {
		corr[i].resize(nFeat);
	}

	features.resize(nSamp * nFeat);
	x.resize(nSamp);
	y.resize(nSamp);

	//get first random samples
	int i = starti;
	Point<T> p;
	while (i < nSamp) {
		clhs.getRandomPoint(p);

		if (points.contains((((uint64_t) p.x) << 32) | ((uint64_t) p.y))) {
			continue;
		}

		x[i] = p.x;
		y[i] = p.y;
		
		for (int f = 0; f < nFeat; f++) {
			T val = p.p_features[f];
			features[(i * nFeat) + f] = val;

			int q = getQuantile<T>(val, quantiles[f]);
			sampleCountPerQuantile[(f * nSamp) + q]++;
			quantilesOfEachSample[(i * nFeat) + f] = q;
		}

		points.insert((((uint64_t) p.x) << 32) | ((uint64_t) p.y));
		i++;
	}

	//define covariance calculation 
	DALHomogenTable table = DALHomogenTable(features.data(), nSamp, nFeat, [](const T *){}, oneapi::dal::data_layout::row_major);
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

	//features of old (before random new index) index
	std::vector<T> oldf(nFeat);

	//begin annealing schedule. If we have a perfect latin hypercube -- or if we pass enough iterations -- stop iterating.
	while (temp > 0 && objQ != 0) {
		size_t i; //the index within the indices, x, y, and features vector so we know what to swap without searching
		if (dist(rng) < 0.5) {
			//50% of the time, choose a random sample to replace
			i = indexDist(rng);
		}
		else {
			//50% of the time, choose the worst sample to replace

			//get the sample with the worst redundancy. In other words, the one which is in the
			//most over-represented quantile across all features
			size_t worstRedundancyIndex = 0;
			size_t worstRedundancy = 0;
			for (size_t si = starti; si < nSamp; si++) {
				size_t curSampleRedundancy = 0;
				for (size_t fi = 0; fi < nFeat; fi++) {
					size_t q = quantilesOfEachSample[(si * nFeat) + fi];
					curSampleRedundancy += sampleCountPerQuantile[(fi * nSamp) + q];
				}
				if (curSampleRedundancy > worstRedundancy) {
					worstRedundancy = curSampleRedundancy;
					worstRedundancyIndex = si;
				}
			}
			i = worstRedundancyIndex;
		}

		//move selected replacement to 'oldf' vector to retain the old values in case we revert back
		//to that state
		std::memcpy(oldf.data(), features.data() + (i * nFeat), nFeat * sizeof(T));

		//select a new index
		Point<T> p;
		clhs.getRandomPoint(p);
		while (points.contains((((uint64_t) p.x) << 32) | ((uint64_t) p.y))) {
			clhs.getRandomPoint(p);
		}

		//move new features into feature vector
		std::memcpy(features.data() + (i * nFeat), p.p_features, nFeat * sizeof(T));
	
		//recalculate sample count per quantile
		std::vector<int> oldq(nFeat);
		std::vector<int> newq(nFeat);
		for (int f = 0; f < nFeat; f++) {
			//decrement based removal of on old features
			int q = getQuantile(oldf[f], quantiles[f]);
			oldq[f] = q;
			sampleCountPerQuantile[(f * nSamp) + q]--;

			//increment based on inputof new features
			q = getQuantile(p.p_features[f], quantiles[f]);
			newq[f] = q;
			sampleCountPerQuantile[(f * nSamp) + q]++;
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

		bool keep = dist(rng) < std::exp(-1 * delta / temp);

		if (keep) {
			//update the new changes
			points.erase((((uint64_t) x[i]) << 32) | ((uint64_t) y[i]));
			points.insert((((uint64_t) p.x) << 32) | ((uint64_t) p.y));

			x[i] = p.x;
			y[i] = p.y;
		
			std::memcpy(quantilesOfEachSample.data() + (i * nFeat), newq.data(), sizeof(int) * nFeat);

			objC = newObjC;
			objQ = newObjQ;
			obj = newObj;
		}
		else {
			//revert back to old changes
			for (int f = 0; f < nFeat; f++) {
				sampleCountPerQuantile[(f * nSamp) + newq[f]]--;
				sampleCountPerQuantile[(f * nSamp) + oldq[f]]++;
			}

			std::memcpy(
				reinterpret_cast<void *>(features.data() + (i * nFeat)),
				reinterpret_cast<void *>(oldf.data()),
				nFeat * sizeof(T)
			);
		}

		//update annealing temperature
		temp -= d;
	}

	//add samples to output layer
	helper::Field fieldExistingFalse("existing", 0);
	for (int i = starti ; i < nSamp; i++) {
		const auto [xCoord, yCoord] = helper::sample_to_point(GT, x[i], y[i]);
		OGRPoint point = OGRPoint(xCoord, yCoord);
		existing.used ? 
			helper::addPoint(&point, p_layer, &fieldExistingFalse) :
			helper::addPoint(&point, p_layer);

		if (plot) {
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
		}
	}
}

/**
 * @ingroup clhs
 * This function conducts Conditioned Latin Hypercube Sampling (CLHS) on a given raster.
 *
 * First, metadata is acquired for each band within the raster, and an output vector
 * file is created in-memory which is where chosen sample points will be written to.
 *
 * The readRaster() function is called, which iterates through the raster in blocks
 * and calculates the quantiles for each feature, as well as the correlation matrix,
 * for the entire sampling space. It also selects pixels and adds them to a pool.
 *
 * The selectSamples() function utilizes a simulated annealing algorithm taking
 * into account the differences in correlation matrix between the sample and the
 * overall sample space, as well as how close to a latin hypercube the sample is.
 * The points which this algorithm choses from are from the pool of pixels added
 * in the readRaster() function.
 *
 * @param GDALRasterWrapper *p_raster
 * @param int nSamp
 * @param int iterations
 * @param GDALVectorWrapper *p_access
 * @param std::string layerName
 * @param double buffInner
 * @param double buffOuter
 * @param bool plot
 * @param std::string tempFolder
 * @param std::string filename
 * @returns std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *>
 */
std::tuple<std::vector<std::vector<double>>, vector::GDALVectorWrapper *>
clhs(
	raster::GDALRasterWrapper *p_raster, 
	int nSamp,
	int iterations,
	vector::GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	vector::GDALVectorWrapper *p_existing,
	size_t replace,
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

	std::vector<helper::RasterBandMetaData> bands(p_raster->getBandCount());
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

	vector::GDALVectorWrapper *p_wrapper = new vector::GDALVectorWrapper(p_samples, std::string(p_raster->getDataset()->GetProjectionRef()));
	OGRLayer *p_layer = p_samples->CreateLayer("samples", p_wrapper->getSRS(), wkbPoint, nullptr);
	if (!p_layer) {
		throw std::runtime_error("unable to create output dataset layer.");
	}

	access::Access access(
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

	existing::Existing existing(
		p_existing,
		p_raster,
		GT,
		width,
		p_layer,
		false,
		xCoords,
		yCoords,
		false	
	);

	//fast random number generator using xoshiro256++
	//https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf
	xso::xoshiro_4x64_plus rng;
	uint64_t multiplier = helper::getProbabilityMultiplier(
		width, 
		height, 
		p_raster->getPixelWidth(), 
		p_raster->getPixelHeight(), 
		8, 
		MILLION * 100, 
		false, 
		access.area
	);
	helper::RandValController rand(bands[0].xBlockSize, bands[0].yBlockSize, multiplier, &rng);

	//get data type for all bands
	GDALDataType type = GDT_Float32;
	for (const helper::RasterBandMetaData& band : bands) {
		if (band.type == GDT_Float64) {
			type = GDT_Float64;
			break;
		}
	}

	if (type == GDT_Float64) {	
		std::vector<std::vector<double>> quantiles;
		
		//create instance of data management class
		CLHSDataManager<double> clhs(nFeat, nSamp, &rng, existing.count());

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRaster<double>(bands, clhs, access, existing, rand, type, quantiles, sizeof(double), width, height, nFeat, nSamp);

		//select samples and add them to output layer
		selectSamples<double>(quantiles, clhs, rng, existing, replace, iterations, nSamp, nFeat, p_layer, GT, plot, xCoords, yCoords);
	}
	else { //type == GDT_Float32	
		std::vector<std::vector<float>> quantiles;

		//create instance of data management class
		CLHSDataManager<float> clhs(nFeat, nSamp, &rng, existing.count());

		//read raster, calculating quantiles, correlation matrix, and adding points to sample from.
		readRaster<float>(bands, clhs, access, existing, rand, type, quantiles, sizeof(float), width, height, nFeat, nSamp);

		//select samples and add them to output layer
		selectSamples<float>(quantiles, clhs, rng, existing, replace, iterations, nSamp, nFeat, p_layer, GT, plot, xCoords, yCoords);
	}

	if (filename != "") {
		try {
			p_wrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	return {{xCoords, yCoords}, p_wrapper};
}

} //namespace clhs
} //namespace sgs
