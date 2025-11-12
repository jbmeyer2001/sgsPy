/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of stratified sampling
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#include <iostream>
#include <random>

#include "access.h"
#include "existing.h"
#include "helper.h"
#include "raster.h"
#include "vector.h"

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <xoshiro.h>

/**
 * Helper function which calculates the count of values for each strata
 * depending on the allocation method.
 *
 * @param std::string allocation method
 * @param std::vector<size_t>& sizes of each strata
 * @param std::vector<double> weights of each strata
 * @param size_t total number of pixels
 * @returns std::vector<size_t> counts of each stratum
 */
std::vector<int64_t>
calculateAllocation(
	int64_t numSamples,
	std::string allocation, 
	std::vector<int64_t> strataCounts, 
	std::vector<double> weights,
	int64_t numPixels)
{
	std::vector<int64_t> retval;
	int64_t remainder = numSamples;
	int64_t numStrata = strataCounts.size();
	if (allocation == "prop") {
		//allocate the samples per stratum according to stratum size
		int64_t pixelsPerSample = numPixels / numSamples;

		//add 1 if pixelsPerSample was truncated down, to avoid adding too many samples
		pixelsPerSample += static_cast<int64_t>(pixelsPerSample * numSamples < numPixels);

		for (int64_t i = 0; i < numStrata; i++) {
			int64_t count = strataCounts[i] / pixelsPerSample;
			retval.push_back(count);
			remainder -= count;
		}
	}
	else if (allocation == "equal") {
		//determine the count of samples per strata
		int64_t strataSampleCount = numSamples / numStrata;
		
		for (int64_t i = 0; i < numStrata; i++) {
			retval.push_back(strataSampleCount);
			remainder -= strataSampleCount;
		}
	}
	else if (allocation == "manual" || allocation == "optim") { 
		//allocate samples accordign to weights. 
		//Optim allocation calculates these weights whereas when using manual the weights are given by the user
		for (int64_t i = 0; i < numStrata; i++) {
			int64_t count = static_cast<int64_t>(static_cast<double>(numSamples) * weights[i]);
			retval.push_back(count);
			remainder -= count;
		}
	}
	else { 
		throw std::runtime_error("allocation method must be one of 'prop', 'equal', 'manual', or 'optim'.");
	}

	//redistribute remainder pixels among strata, and check strata sizes
	for (int64_t i = numStrata; i > 0; i--) {
		int64_t extra = remainder / i;
		retval[i - 1] += extra;
		remainder -= extra;

		if (retval[i - 1] > strataCounts[i - 1]) {
			std::cout << "warning: strata " << i - 1 << " does not have enough pixels for the full " << retval[i - 1] << " samples it should recieve. There will be less than " << numSamples << " final samples." << std::endl;
			retval[i - 1] = strataCounts[i - 1];
		}
	}

	return retval;
}

/**
 *
 */
struct OptimAllocationDataManager {
	RasterBandMetaData band;
	std::vector<Variance> variances;
	bool used = false;

	/**
	 *
	 */
	OptimAllocationDataManager(std::string allocation) {
		this->used = allocation == "optim";
	}

	/**
	 * Copy constructor.
	 */
	OptimAllocationDataManager(const OptimAllocationDataManager &other) {
		this->band = other.band;
		this->variances = other.variances;
		this->used = other.used;
	}

	/**
	 *
	 */
	~OptimAllocationDataManager() {
		if (this->used) {
			VSIFree(this->band.p_buffer);
		}
	}

	/**
	 *
	 */
	inline void
	init(RasterBandMetaData& band, int numStrata) {
		this->band = band;
		this->variances.resize(numStrata);
	}

	/**
	 *
	 */
	inline void
	readNewBlock(int xBlockSize, int yBlockSize, int xBlock, int yBlock, int xValid, int yValid, void *p_buffer) {
		rasterBandIO(this->band, p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true);
	}

	/**
	 *
	 */
	inline void
	update(int index, int strata, void *p_buffer) {
		double val = getPixelValueDependingOnType<double>(this->band.type, p_buffer, index);
		variances[strata].update(val);	
	}

	/**
	 *
	 */
	inline void
	update(const OptimAllocationDataManager &other) {
		for (size_t i = 0; i < variances.size(); i++) {
			this->variances[i].update(other.variances[i]);
		}
	}

	/**
	 *
	 */
	inline std::vector<double>
	getAllocationPercentages(void) {
		std::vector<double> retval(variances.size());

		double total = 0;
		for (size_t i = 0; i < variances.size(); i++) {
			double stdev = variances[i].getStdev();
			double count = static_cast<double>(variances[i].getCount());
			double product = count == 0 ? 0 : stdev * count; //if count == 0 stdev will be nan, this stops the nan from spreading
			retval[i] = product;
			total += product;
		}

		for (size_t i = 0; i < variances.size(); i++) {
			retval[i] = retval[i] / total;
		}

		return retval;
	}
};

inline uint64_t
getProbabilityMultiplier(GDALRasterWrapper *p_raster, int numSamples, bool useMindist, double accessibleArea, bool queinnec) {
	double height = static_cast<double>(p_raster->getHeight());
	double width = static_cast<double>(p_raster->getWidth());
	double samples = static_cast<double>(numSamples);

	
	double numer = !queinnec ? 
		samples * 4 * (useMindist ? 3 : 1) :
		samples * 32 * (useMindist ? 3 : 1); //much more likely to add a queinnec pixel
	double denom = height * width;

	if (accessibleArea != -1) {
		double pixelHeight = static_cast<double>(p_raster->getPixelHeight());
		double pixelWidth = static_cast<double>(p_raster->getPixelWidth());
		double totalArea = width * pixelWidth * height * pixelHeight;

		numer *= (totalArea / accessibleArea);
	}

	uint8_t bits = static_cast<uint8_t>(std::ceil(std::log2(denom) - std::log2(numer)));
	return (bits <= 0) ? 0 : (1 << bits) - 1;
}

/**
 *
 */
class IndexStorageVectors {
private:
	std::vector<int64_t> strataCounts;

	std::vector<int64_t> indexCountPerStrata;
	std::vector<std::vector<Index>> indexesPerStrata;
	
	std::vector<int64_t> firstXIndexCountPerStrata;
	std::vector<std::vector<Index>> firstXIndexesPerStrata;

	int64_t numStrata;
	int64_t x;

public:
	/**
	 *
	 */
	IndexStorageVectors(int64_t numStrata, int64_t x) {
		this->numStrata = numStrata;
		this->x = x;

		this->strataCounts.resize(numStrata);
		
		this->indexCountPerStrata.resize(numStrata, 0);
		this->indexesPerStrata.resize(numStrata);
		
		this->firstXIndexCountPerStrata.resize(numStrata, 0);
		this->firstXIndexesPerStrata.resize(numStrata);
		for (int64_t i = 0; i < numStrata; i++) {
			this->firstXIndexesPerStrata[i].resize(x);
		}
	}

	/**
	 * Copy constructor
	 */
	IndexStorageVectors(const IndexStorageVectors &other) {
		this->strataCounts = other.strataCounts;
		this->indexCountPerStrata = other.indexCountPerStrata;
		this->indexesPerStrata = other.indexesPerStrata;
		this->firstXIndexCountPerStrata = other.firstXIndexCountPerStrata;
		this->firstXIndexesPerStrata = other.firstXIndexesPerStrata;
		this->numStrata = other.numStrata;
		this->x = other.x;
	}

	/**
	 *
	 */
	inline void
	update(const IndexStorageVectors &other) {
		for (int64_t i = 0; i < numStrata; i++) {
			this->strataCounts[i] += other.strataCounts[i];

			if (other.indexCountPerStrata[i] > 0) {
				size_t thisSize = this->indexesPerStrata[i].size();
				size_t otherSize = other.indexesPerStrata[i].size();
				size_t totalSize = thisSize + otherSize;
				this->indexesPerStrata[i].resize(totalSize);
				
				//copy values from other indices to this indices
				for (size_t j = 0; j < otherSize; j++) {
					this->indexesPerStrata[i][thisSize + j] = other.indexesPerStrata[i][j];
				}
				
				this->indexCountPerStrata[i] += other.indexCountPerStrata[i];
			}

			if (other.firstXIndexCountPerStrata[i] > 0 && this->firstXIndexCountPerStrata[i] <= this->x) {
				int64_t thisCount = this->firstXIndexCountPerStrata[i];
				int64_t otherCount = other.firstXIndexCountPerStrata[i];
				int64_t total = thisCount + otherCount;
				if (total > this->x) {
					std::vector<Index>().swap(this->firstXIndexesPerStrata[i]);
					this->firstXIndexCountPerStrata[i] = this->x + 1;
				}
				else {
					//copy values from other indices to this indices
					for (int64_t j = 0; j < otherCount; j++) {
						this->firstXIndexesPerStrata[i][thisCount + j] = other.firstXIndexesPerStrata[i][j];
					}

					this->firstXIndexCountPerStrata[i] = total;
				}
			}
		}
	}


	/**
	 *
	 */
	inline void
	updateStrataCounts(int val) {
		this->strataCounts[val]++;
	}

	/**
	 *
	 */
	inline void
	updateIndexesVector(int val, Index& index) {
		this->indexesPerStrata[val].push_back(index);
		this->indexCountPerStrata[val]++;
	}

	/**
	 *
	 */
	inline void
	updateFirstXIndexesVector(int val, Index& index) {
		int64_t i = firstXIndexCountPerStrata[val];
		if (i < this->x) {
			firstXIndexesPerStrata[val][i] = index;
			firstXIndexCountPerStrata[val]++;	
		}
		else if (i == this->x) {
			std::vector<Index>().swap(firstXIndexesPerStrata[val]);
			firstXIndexCountPerStrata[val]++;
		}
	}
	
	/**
	 *
	 */
	inline std::vector<std::vector<Index> *>
	getStrataIndexVectors(std::vector<int64_t> existing, std::vector<int64_t> strataSampleCounts, xso::xoshiro_4x64_plus& rng) {
		std::vector<std::vector<Index> *> retval(numStrata);

		for (int64_t i = 0; i < numStrata; i++) {
			int64_t existingSamples = existing[i];
			int64_t desiredSamples = strataSampleCounts[i];
			int64_t remainingSamples = desiredSamples - existingSamples;

			int64_t probIndexesCount = indexCountPerStrata[i];
			int64_t firstXIndexesCount = firstXIndexCountPerStrata[i];

			if (probIndexesCount >= remainingSamples || firstXIndexesCount > x) {
				auto begin = this->indexesPerStrata[i].begin();
				auto end = this->indexesPerStrata[i].end();
				std::shuffle(begin, end, rng);
				retval[i] = &this->indexesPerStrata[i];
			}
			else {
				this->firstXIndexesPerStrata[i].resize(firstXIndexesCount);
				auto begin = this->firstXIndexesPerStrata[i].begin();
				auto end = this->firstXIndexesPerStrata[i].end();
				std::shuffle(begin, end, rng);
				retval[i] = &this->firstXIndexesPerStrata[i];
			}
		}

		return retval;
	}

	/**
	 *
	 */
	inline std::vector<int64_t>
	getStrataCounts(void) {
		return this->strataCounts;
	}

	/**
	 *
	 */
	inline int64_t
	getNumDataPixels(void) {
		int64_t retval = 0;

		for (const int64_t& count : this->strataCounts) {
			retval += count;
		}

		return retval;
	}
};

/**
 *
 */
class RandValController {
private:
	std::vector<bool> randVals;
	size_t randValIndex = 0;
	uint64_t multiplier = 0;
	xso::xoshiro_4x64_plus *p_rng = nullptr;

public:
	/**
	 *
	 */
	RandValController(int xBlockSize, int yBlockSize, uint64_t multiplier, xso::xoshiro_4x64_plus *p_rng) {
		this->randVals.resize(xBlockSize * yBlockSize);
		this->randValIndex = static_cast<size_t>(xBlockSize * yBlockSize);
		this->multiplier = multiplier;
		this->p_rng = p_rng;
	}

	/**
	 *
	 */
	inline void 
	calculateRandValues(void) {
		for (size_t i = 0; i < randValIndex; i++) {
			randVals[i] = (((*p_rng)() >> 11) & multiplier) == multiplier;
		}
		randValIndex = 0;
	}

	/**
	 *
	 */
	inline bool 
	next(void) {
		bool retval = randVals[randValIndex];
		randValIndex++;
		return retval;
	}
};

/**
 *
 */
struct Mutexes {
	std::mutex band;
	std::mutex rng;
	std::mutex access;
	std::mutex optim;
	std::mutex indices;
};

/**
 *
 */
template <typename T>
std::vector<int64_t>
processBlocksStratRandom(
	int numSamples,
	int numStrata,
	RasterBandMetaData& band,
	Access& access,
	Existing& existing,
	IndexStorageVectors& indices,
	std::vector<std::vector<OGRPoint>>& existingSamples,
	uint64_t multiplier,
	xso::xoshiro_4x64_plus& rng,
	std::string allocation,
	OptimAllocationDataManager& optim,
	std::vector<double> weights,
	int width,
	int height,
	Mutexes& mutexes,
	int threads) 
{
	T nanInt = static_cast<T>(band.nan);

	int xBlockSize = band.xBlockSize;
       	int yBlockSize = band.yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	//TODO CHANGE THIS BACK
	//int chunkSize = yBlocks / threads;
	int chunkSize = yBlocks / 8;

	if (chunkSize == 0) {
		chunkSize = yBlocks;
	}
	
	//create thread pool and acquire Python's global interpreter lock (GIL)
	pybind11::gil_scoped_acquire acquire;
	boost::asio::thread_pool pool(threads);

	for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
		int yBlockEnd = std::min(yBlockStart + chunkSize, yBlocks);

		boost::asio::post(pool, [
			xBlockSize,
			yBlockSize,
			yBlockStart,
			yBlockEnd,
			xBlocks,
			multiplier,
			nanInt,
			&band,
			&access,
			&existing,
			&existingSamples,
			&indices,
			&rng,
			&optim,
			&mutexes
		] {
			//create rand val controller for this thread
			RandValController rand(xBlockSize, yBlockSize, multiplier, &rng);
	
			//allocate buffers 
			T *p_buffer = reinterpret_cast<T *>(VSIMalloc3(xBlockSize, yBlockSize, band.size));
			int8_t *p_access = nullptr;
			void *p_optim = nullptr;
			if (access.used) {
				p_access = reinterpret_cast<int8_t *>(VSIMalloc3(xBlockSize, yBlockSize, access.band.size));
			}
			if (optim.used) {
				p_optim = VSIMalloc3(xBlockSize, yBlockSize, optim.band.size);
			}	

			//create thread-specific objects, so that
			//we can update the original at the end of this threads execution
			//rather than acquiring a lock every time we want to edit one
			//of these objects during execution, which would have significant overhead.
			OptimAllocationDataManager threadOptim(optim);
			IndexStorageVectors threadIndices(indices);
			std::vector<std::vector<OGRPoint>> threadExistingSamples(existingSamples.size());

			for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
				for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
					int xValid, yValid;
		
					//read block
					band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
					rasterBandIO(band, reinterpret_cast<void *>(p_buffer), xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true);

					//read access block
					if (access.used) {
						rasterBandIO(access.band, reinterpret_cast<void *>(p_access), xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true);
					}

					//read mraster block for optim
					if (threadOptim.used) {
						threadOptim.readNewBlock(xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, p_optim);
					}

					//calculate rand vals
					mutexes.rng.lock();
					rand.calculateRandValues();
					mutexes.rng.unlock();

					//iterate through block and update vectors
					for (int y = 0; y < yValid; y++) {
						int blockIndex = y * xBlockSize;
						for (int x = 0; x < xValid; x++) {
							T val = p_buffer[blockIndex];
							Index index = {x + xBlock * xBlockSize, y + yBlock * yBlockSize};

							bool isNan = val == nanInt;
							bool accessible = !access.used || p_access[blockIndex] != 1;
							bool alreadySampled = existing.used && existing.containsIndex(index.x, index.y);

							//check nan
							if (isNan) {
								blockIndex++;
								continue;
							}

							//update optim allocation variance calculations
							if (threadOptim.used) {
								threadOptim.update(blockIndex, val, p_optim);
							}

							//udpate strata counts
							threadIndices.updateStrataCounts(val);
					
							//update existing sampled strata
							if (alreadySampled) {
								threadExistingSamples[val].push_back(existing.getPoint(index.x, index.y));
							}

							//add val to stored indices
							if (accessible && !alreadySampled) {
								threadIndices.updateFirstXIndexesVector(val, index);

								if (rand.next()) {
									threadIndices.updateIndexesVector(val, index);
								}					
							}

							//increment block index
							blockIndex++;
						}
					}
				}
			}

			//free memory allocated by this thread
			VSIFree(p_buffer);
			if (access.used) {
				VSIFree(p_access);
			}
			if (optim.used) {
				VSIFree(p_optim);
			}

			//update overall optim object with thread-specific optim object
			mutexes.optim.lock();
			optim.update(threadOptim);
			mutexes.optim.unlock();

			//update overall indices object with thread-specific indices object
			mutexes.indices.lock();
			indices.update(threadIndices);
			mutexes.indices.unlock();

			//update existing samples with thread-specific vector
			if (existing.used && threadExistingSamples.size() > 0) {
				mutexes.indices.lock();
				existingSamples.resize(existingSamples.size() + threadExistingSamples.size());

				std::memcpy(
					(void *)((size_t)existingSamples.data() + existingSamples.size() * sizeof(Index)), //dest
					(void *)threadExistingSamples.data(),						   //src
					threadExistingSamples.size() * sizeof(Index)					   //count
				);
				mutexes.indices.unlock();
			}
		});
	}
	
	//join thread pool and release Python's global interpreter lock (GIL)
	pool.join();
	pybind11::gil_scoped_release release;

	if (optim.used) {
		weights = optim.getAllocationPercentages();
	}

	std::vector<int64_t> strataSampleCounts = calculateAllocation(
		numSamples,
		allocation,
		indices.getStrataCounts(),
		weights,
		indices.getNumDataPixels()
	);

	return strataSampleCounts;
}

/**
 *
 */
struct FocalWindow {
	int wrow, wcol;
	int width;
	int vpad, hpad;
	std::vector<bool> m;
	std::vector<bool> valid;

	FocalWindow(int wrow, int wcol, int width) {
		this->wrow = wrow;
		this->vpad = wrow / 2;
		this->wcol = wcol;
		this->hpad = wcol / 2;
		this->width = width;
		this->m.resize(wrow * width, false);
		this->valid.resize(wrow * width, false);
	}

	inline void
	reset(int row) {
		int start = (row % wrow) * this->width;
		int end = start + this->width;
		for (int i = start; i < end; i++) {
			m[i] = false;
			valid[i] = false;
		}
	}

	inline bool
	check(int x, int y) {
		y = y % wrow;
		switch (wrow) {
			case 3:
				return this->m[x] && 
				       this->m[x + width] && 
				       this->m[x + width * 2] &&
				       this->valid[x + y * width];		
			case 5:
				return this->m[x] && 
				       this->m[x + width] && 
				       this->m[x + width * 2] && 
				       this->m[x + width * 3] && 
				       this->m[x + width * 4] &&
				       this->valid[x + y * width];
			case 7:
				return this->m[x] && 
				       this->m[x + width] && 
				       this->m[x + width * 2] && 
				       this->m[x + width * 3] && 
				       this->m[x + width * 4] &&
				       this->m[x + width * 5] &&
				       this->m[x + width * 6] && 
				       this->valid[x + y * width];
			default:
				throw std::runtime_error("wrow must be one of 3, 5, 7.");				
		}
	}
};

/**
 *
 */
template <typename T>
std::vector<int64_t>
processBlocksStratQueinnec(
	int numSamples,
	int numStrata,
	RasterBandMetaData &band,
	Access& access,
	Existing& existing,
	IndexStorageVectors& indices,
	IndexStorageVectors& queinnecIndices,
	std::vector<std::vector<OGRPoint>>& existingSamples,
	uint64_t multiplier,
	uint64_t queinnecMultiplier,
	xso::xoshiro_4x64_plus& rng,
	std::string allocation,
	OptimAllocationDataManager& optim,
	std::vector<double> weights,
	int width,
	int height,
	int wrow,
	int wcol,
	Mutexes& mutexes,
	int threads) 
{
	T nanInt = static_cast<T>(band.nan);

	//adjust blocks to be a large chunk of scanlines
	int xBlockSize = band.xBlockSize;
	int yBlockSize = band.yBlockSize;
	if (xBlockSize != width) {
		xBlockSize = width;
	}
	else {
		yBlockSize = 128;
	}

	int yBlocks = (height + yBlockSize - 1) / yBlockSize;
	int chunkSize = yBlocks / threads;

	if (chunkSize == 0) {
		chunkSize = yBlocks;
	}

	//create thread pool and acquire Python's global interpreter lock (GIL)
	pybind11::gil_scoped_acquire acquire;
	boost::asio::thread_pool pool(threads);

	for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
		int yBlockEnd = std::min(yBlockStart + chunkSize, yBlocks);

		boost::asio::post(pool, [
			wrow,
			wcol,
			height,
			width,
			xBlockSize,
			yBlockSize,
			yBlockStart,
			yBlockEnd,
			multiplier,
			queinnecMultiplier,
			nanInt,
			&band,
			&access,
			&existing,
			&indices,
			&queinnecIndices,
			&existingSamples,
			&rng,
			&optim,
			&mutexes
		] {
			//create rand val controllers for this thread
			RandValController rand(xBlockSize, yBlockSize, multiplier, &rng);
			RandValController queinnecRand(xBlockSize, yBlockSize, queinnecMultiplier, &rng);	

			//allocate buffers
			T *p_buffer = reinterpret_cast<T *>(VSIMalloc3(xBlockSize, yBlockSize, band.size));
		       	int8_t *p_access = nullptr;
			void *p_optim = nullptr;
			if (access.used) {
				p_access = reinterpret_cast<int8_t *>(VSIMalloc3(xBlockSize, yBlockSize, access.band.size));
			}	
			if (optim.used) {
				p_optim = VSIMalloc3(xBlockSize, yBlockSize, optim.band.size);
			}
			
			//create thread-specific objects, so that
			//we can update the original at the end of this threads execution
			//rather than acquiring a lock every time we want to edit one
			//of these objects during execution, which would have significant overhead.
			OptimAllocationDataManager threadOptim(optim);
			IndexStorageVectors threadIndices(indices);
			std::vector<std::vector<OGRPoint>> threadExistingSamples(existingSamples.size());

			FocalWindow fw(wrow, wcol, width);

			int fwyi = 0;
			for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
				//read block
				int xOff = 0;
				int yOff = yBlock * yBlockSize;
				int xValid = width;
				int yValid = std::min(yBlockSize, height - yBlock * yBlockSize);
		
				//read block
				mutexes.band.lock();
				if (yBlock == 0) { 
					band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, band.p_buffer, xValid, yValid, band.type, 0, 0);
				}
				else {
					int stratYOff = yOff - fw.vpad * 2;
					int stratYValid = yValid + fw.vpad * 2;
					band.p_band->RasterIO(GF_Read, xOff, stratYOff, xValid, stratYValid, reinterpret_cast<void *>(p_buffer), xValid, stratYValid, band.type, 0, 0);
				}
				mutexes.band.unlock();
		
				//read access block
				if (access.used) {
					mutexes.access.lock();
					access.band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, reinterpret_cast<void *>(p_access), xValid, yValid, GDT_Int8, 0, 0);
					mutexes.access.unlock();
				}
		
				//read mraster block for optim
				if (threadOptim.used) {
					mutexes.optim.lock();
					threadOptim.band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, p_optim, xValid, yValid, optim.band.type, 0, 0);
					mutexes.optim.unlock();
				}
			
				//calculate rand vals
				mutexes.rng.lock();
				rand.calculateRandValues();
				queinnecRand.calculateRandValues();
				mutexes.rng.unlock();
		
				//calculate the within-block index. In the case where there is a padding at the top of the block,
				//adjust for that.
				int newBlockStart = (yBlock == 0) ? 0 : width * fw.vpad * 2;
		
				//iterate through block and update vectors
				for (int y = 0; y < yValid; y++) {
					//reset upcomming row in focal window matrix
					fw.reset(yBlock * yBlockSize + y);	
					
					//the indexes from 0 to the horizontal pad cannot be a queinnec index
					//because they are too close to the edges of the raster	
					for (int x = 0; x < fw.hpad; x++) {
						T val = p_buffer[newBlockStart + y * width + x];
						Index index = {x, y + yBlock * yBlockSize};
		
						bool isNan = val == nanInt;
						bool accessible = !access.used || p_access[y * width + x] != 1;
						bool alreadySampled = existing.used && existing.containsIndex(index.x, index.y);
		
						//check nan
						if (isNan) {
							continue;
						}
		
						//update optim allocation variance calculations
						if (threadOptim.used) {
							threadOptim.update(y * width + x, val, p_optim);
						}
		
						//update strata counts
						threadIndices.updateStrataCounts(val);
		
						//update existing samled strata
						if (alreadySampled) {
							existingSamples[val].push_back(existing.getPoint(index.x, index.y));
						}
		
						//add val to stored indices
						if (accessible && !alreadySampled) {
							threadIndices.updateFirstXIndexesVector(val, index);
		
							if (rand.next()) {
								threadIndices.updateIndexesVector(val, index);
							}
						}
					}
			
					for (int x = fw.hpad; x < width - fw.hpad; x++) {
						T val = p_buffer[newBlockStart + y * width + x];
						Index index = {x, y + yBlock * yBlockSize};
		
						bool isNan = val == nanInt;
						bool accessible = !access.used || p_access[y * width + x] != 1;
						bool alreadySampled = existing.used && existing.containsIndex(index.x, index.y);
						
						//check nan
						if (isNan) {
							continue;
						}
		
						//update optim allocation variance calculations
						if (threadOptim.used) {
							threadOptim.update(y * width + x, val, p_optim);
						}
		
						//update strata counts
						threadIndices.updateStrataCounts(val);
		
						//update exisitng sampled strata
						if (alreadySampled) {
							existingSamples[val].push_back(existing.getPoint(index.x, index.y));
						};
		
						//set focal window matrix value by checking horizontal indices within focal window
						int start = newBlockStart + y * width + x - fw.hpad;
						switch (fw.wcol) {
							case 3:
								fw.m[fwyi + x] = p_buffer[start] == p_buffer[start + 1] &&
										p_buffer[start] == p_buffer[start + 2];
								break;
							case 5: 
								fw.m[fwyi + x] = p_buffer[start] == p_buffer[start + 1] &&
										p_buffer[start] == p_buffer[start + 2] &&
										p_buffer[start] == p_buffer[start + 3] &&
										p_buffer[start] == p_buffer[start + 4];
								break;
							case 7: 
								fw.m[fwyi + x] = p_buffer[start] == p_buffer[start + 1] &&
										p_buffer[start] == p_buffer[start + 2] &&
										p_buffer[start] == p_buffer[start + 3] &&
										p_buffer[start] == p_buffer[start + 4] &&
										p_buffer[start] == p_buffer[start + 5] &&
										p_buffer[start] == p_buffer[start + 6];
								break;
							default:
								throw std::runtime_error("wcol must be one of 3, 5, 6.");
						}
						
						//add val to stored indices and update focal window validity	
						if (accessible && !alreadySampled) {
							threadIndices.updateFirstXIndexesVector(val, index);
		
							if (rand.next()) {
								threadIndices.updateIndexesVector(val, index);
							}
							
							fw.valid[fwyi + x] = true;
						}
		
						//add index to queinnec indices if surrounding focal window is all the same
						Index fwIndex = {x, index.y - fw.vpad};
						if (fw.check(fwIndex.x, fwIndex.y)) {
							bool add = false;
							int start = newBlockStart + y * width + x - 2 * width * fw.vpad;
							switch (fw.wrow) {
								case 3:
									add = p_buffer[start] == p_buffer[start + width * 1] &&
										p_buffer[start] == p_buffer[start + width * 2];
									break;
								case 5:
									add = p_buffer[start] == p_buffer[start + width * 1] &&
										p_buffer[start] == p_buffer[start + width * 2] &&
										p_buffer[start] == p_buffer[start + width * 3] &&
										p_buffer[start] == p_buffer[start + width * 4];
									break;
								case 7:
									add = p_buffer[start] == p_buffer[start + width * 1] &&
										p_buffer[start] == p_buffer[start + width * 2] &&
										p_buffer[start] == p_buffer[start + width * 3] &&
										p_buffer[start] == p_buffer[start + width * 4] &&
										p_buffer[start] == p_buffer[start + width * 5] &&
										p_buffer[start] == p_buffer[start + width * 6];
									break;
							}
		
							if (add) {
								//(we know if we've made it here that val is the same for both the fw add and the current index)
								queinnecIndices.updateFirstXIndexesVector(val, fwIndex);
		
								if (queinnecRand.next()) {
									queinnecIndices.updateIndexesVector(val, fwIndex);
								}
							}
						}
					}
		
					//the indexes from width - horizontal pad to width cannot be a queinnec index
					//because they are too close to the edges of the raster	
					for (int x = width - fw.hpad; x < width; x++) {
						T val = p_buffer[newBlockStart + y * width + x];
						Index index = {x, y + yBlock * yBlockSize};
		
						bool isNan = val == nanInt;
						bool accessible = !access.used || p_access[y * width + x] != 1;
						bool alreadySampled = existing.used && existing.containsIndex(index.x, index.y);
		
						//check nan
						if (isNan) {
							continue;
						}
		
						//update optim allocation variance calculations
						if (threadOptim.used) {
							threadOptim.update(y * width + x, val, p_optim);
						}
		
						//update strata counts
						threadIndices.updateStrataCounts(val);
		
						//update existing sampled strata
						if (alreadySampled) {
							existingSamples[val].push_back(existing.getPoint(index.x, index.y));
						}
		
						//add val to stored indices
						if (accessible && !alreadySampled) {
							threadIndices.updateFirstXIndexesVector(val, index);
		
							if (rand.next()) {
								threadIndices.updateIndexesVector(val, index);
							}
						}
					}
		
					fwyi += width;
					if (fwyi == width * fw.wrow) {
						fwyi = 0;
					}
				}
			}

			
			//free memory allocated by this thread
			VSIFree(p_buffer);
			if (access.used) {
				VSIFree(p_access);
			}
			if (optim.used) {
				VSIFree(p_optim);
			}

			//update overall optim object with thread-specific optim object
			mutexes.optim.lock();
			optim.update(threadOptim);
			mutexes.optim.unlock();

			//update overall indices object with thread-specific indices object
			mutexes.indices.lock();
			indices.update(threadIndices);
			mutexes.indices.unlock();

			//update existing samples with thread-specific vector
			if (existing.used && threadExistingSamples.size() > 0) {
				mutexes.indices.lock();
				existingSamples.resize(existingSamples.size() + threadExistingSamples.size());

				std::memcpy(
					(void *)((size_t)existingSamples.data() + existingSamples.size() * sizeof(Index)), //dest
					(void *)threadExistingSamples.data(),						   //src
					threadExistingSamples.size() * sizeof(Index)					   //count
				);
				mutexes.indices.unlock();
			}
	
		});
	}	
	
	//join thread pool and release Python's global interpreter lock (GIL)
	pool.join();
	pybind11::gil_scoped_release release;

	if (optim.used) {
		weights = optim.getAllocationPercentages();
	}
	std::vector<int64_t> strataSampleCounts = calculateAllocation(
		numSamples,
		allocation,
		indices.getStrataCounts(),
		weights,
		indices.getNumDataPixels()
	);

	return strataSampleCounts;	
}

/**
 * This function conducts stratified random sampling on the provided stratified raster.
 *
 *
 * First, the raster is iterated over, the no data pixels, and the pixels which
 * are inaccessible are ignored. The remaining accessable data pixel indexes
 * are placed into vectors corresponding to their strata. After this iteration,
 * there will be a a number of vectors equivalent to the number of strata,
 * correspondin strata.
 *
 *
 * Nextl the calculate_allocation function is used to determine the the total
 * number of samples which should be allocated to each strata, depending on the
 * number of pixels in each strata and the allocation method specified by the
 * user.
 *
 *
 * Finally, random indexes in the strata vectors are selected, and using
 * their index value in the stratified raster and geotransform, a point
 * geometry is calculated.
 *
 * When a stratum is allocated more than half the number of total pixels the 
 * strata vector contains, random indexes in the strat vectors are selected 
 * which WONT be included, and the remaining are iterated over in a random order.
 *
 * When mindist is inequal to zero, three times the number of samples are 
 * randomly selected in order to ensure a number of samples as close to the 
 * desired number are selected, as having a large mindist may mean some samples
 * can't be included due to their proximity to other already-selected pixels. 
 */
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
strat(
	GDALRasterWrapper *p_raster,
	int bandNum,
	int64_t numSamples,
	int64_t numStrata,
	std::string allocation,
	std::vector<double> weights,
	GDALRasterWrapper *p_mraster,
	int mrastBandNum,
	std::string method,
	int wrow,
	int wcol,
	double mindist,
	GDALVectorWrapper *p_existing,
	bool force,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string filename,
	std::string tempFolder,
	int threads)
{
	GDALAllRegister();

	bool useMindist = mindist != 0;
	int width = p_raster->getWidth();
	int height = p_raster->getHeight();
	double *GT = p_raster->getGeotransform();

	Mutexes mutexes;

	//step 1: get raster band
	RasterBandMetaData band = p_raster->getRasterBandMetaData(bandNum);
	band.p_mutex = &mutexes.band;
	
	printTypeWarningsForInt32Conversion(band.type);
	
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
		band.xBlockSize,
		band.yBlockSize
	);
	access.band.p_mutex = &mutexes.access;

	std::vector<double> xCoords, yCoords;
	std::vector<std::vector<OGRPoint>> existingSamples(numStrata);	
	Existing existing(
		p_existing,
		GT,
		width,
		nullptr,
		false,
		xCoords,
		yCoords	
	);

	//fast random number generator using xoshiro256+
	//https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf
	xso::xoshiro_4x64_plus rng; 
	uint64_t multiplier = getProbabilityMultiplier(p_raster, numSamples, useMindist, access.area, false);
	uint64_t queinnecMultiplier = getProbabilityMultiplier(p_raster, numSamples, useMindist, access.area, true);

	//normal indices used with both methods, queinnec indices used only with queinnec method
	IndexStorageVectors indices(numStrata, 10000);
	IndexStorageVectors queinnecIndices(numStrata, 10000);	

	OptimAllocationDataManager optim(allocation);
	if (optim.used) {
		RasterBandMetaData optimBand = p_raster->getRasterBandMetaData(bandNum);
		optimBand.p_mutex = &mutexes.optim;
		optim.init(optimBand, numStrata);
	}

	std::vector<int64_t> strataSampleCounts; 
	if (method == "random") {
		switch (band.type) {
			case GDT_Int8:
				strataSampleCounts = processBlocksStratRandom<int8_t>(numSamples, numStrata, band, access, 
										      existing, indices, existingSamples,
							 			      multiplier, rng, allocation, optim,
										      weights, width, height, mutexes, threads);
				break;
			case GDT_Int16:
				strataSampleCounts = processBlocksStratRandom<int16_t>(numSamples, numStrata, band, access, 
										      existing, indices, existingSamples,
							 			      multiplier, rng, allocation, optim,
										      weights, width, height, mutexes, threads);
				break;
			default:
				strataSampleCounts = processBlocksStratRandom<int32_t>(numSamples, numStrata, band, access, 
										      existing, indices, existingSamples,
							 			      multiplier, rng, allocation, optim,
										      weights, width, height, mutexes, threads);
				break;
		}
	}
	else { //method == queinnec
		switch (band.type) {
			case GDT_Int8:
				strataSampleCounts = processBlocksStratQueinnec<int8_t>(numSamples, numStrata, band, access, existing, 
											indices, queinnecIndices, existingSamples, multiplier, 
											queinnecMultiplier, rng, allocation, optim, weights, 
											width, height, wrow, wcol, mutexes, threads);
				break;
			case GDT_Int16:
				strataSampleCounts = processBlocksStratQueinnec<int16_t>(numSamples, numStrata, band, access, existing, 
											indices, queinnecIndices, existingSamples, multiplier, 
											queinnecMultiplier, rng, allocation, optim, weights, 
											width, height, wrow, wcol, mutexes, threads);
				break;
			default:
				strataSampleCounts = processBlocksStratQueinnec<int32_t>(numSamples, numStrata, band, access, existing, 
											indices, queinnecIndices, existingSamples, multiplier, 
											queinnecMultiplier, rng, allocation, optim, weights, 
											width, height, wrow, wcol, mutexes, threads);
				break;
		}
	}

	std::vector<std::vector<Index> *> strataIndexVectors;
	std::vector<size_t> nextIndexes(numStrata, 0);

	std::vector<bool> completedStrata(numStrata, false);
	std::vector<bool> completedStrataQueinnec(numStrata, false);
	std::vector<int64_t> samplesAddedPerStrata(numStrata, 0);
 	int64_t numCompletedStrata = 0;
	int64_t numCompletedStrataQueinnec = 0;
	int64_t curStrata = 0;
	int64_t addedSamples = 0;

	//add existing sample plots
	if (existing.used) {
		if (force) {
			//if force is used, add all samples no matter what
			for (size_t i = 0; i < existingSamples.size(); i++) {
				std::vector<OGRPoint> samples = existingSamples[i];
			       	for (const OGRPoint& point : samples) {
					addPoint(&point, p_layer);

					addedSamples++;
					samplesAddedPerStrata[i]++;

					if (plot) {
						xCoords.push_back(point.getX());
						yCoords.push_back(point.getY());
					}
				}	

				if (samplesAddedPerStrata[i] >= strataSampleCounts[i]) {
					numCompletedStrata++;
					numCompletedStrataQueinnec++;
					completedStrata[i] = true;
					completedStrataQueinnec[i] = true;
				}
			}
		}
		else { //force == false
			for (size_t i = 0; i < existingSamples.size(); i++) {
				std::vector<OGRPoint> samples = existingSamples[i];
				std::shuffle(samples.begin(), samples.end(), rng);
				
				size_t j = 0;
				while (samplesAddedPerStrata[i] < strataSampleCounts[i] && j < samples.size()) {
					OGRPoint point = samples[j];

					if (mindist != 0 && p_layer->GetFeatureCount() != 0) {
						bool add = true;
						for (const auto &p_feature : *p_layer) {
							OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
							if (point.Distance(p_point) < mindist) {
								add = false;
								break;
							}
						}

						if (!add) {
							j++;
							continue;
						}
					}

					addPoint(&point, p_layer);
					
					addedSamples++;
					samplesAddedPerStrata[i]++;
					j++;

					if (plot) {
						xCoords.push_back(point.getX());
						yCoords.push_back(point.getY());
					}
				}

				if (samplesAddedPerStrata[i] == strataSampleCounts[i]) {
					numCompletedStrata++;
					numCompletedStrataQueinnec++;
					completedStrata[i] = true;
					completedStrataQueinnec[i] = true;
				}
			}
		}
	}	

	if (method == "Queinnec") {
		strataIndexVectors = queinnecIndices.getStrataIndexVectors(samplesAddedPerStrata, strataSampleCounts, rng);	

		int64_t curStrata = 0;
		while (numCompletedStrataQueinnec < numStrata && addedSamples < numSamples) {
			if (curStrata == numStrata) {
				curStrata = 0;
			}
			if (completedStrataQueinnec[curStrata]) {
				curStrata++;
				continue;
			}

			int64_t sampleCount = strataSampleCounts[curStrata];
			int64_t samplesAdded = samplesAddedPerStrata[curStrata];
			if (samplesAdded == sampleCount) {
				numCompletedStrataQueinnec++;
				completedStrataQueinnec[curStrata] = true;
				numCompletedStrata++;
				completedStrata[curStrata] = true;
				curStrata++;
				continue;
			}

			std::vector<Index> *strataIndexes = strataIndexVectors[curStrata];
			size_t nextIndex = nextIndexes[curStrata];	
			if (strataIndexes->size() == nextIndex) {
				numCompletedStrataQueinnec++;
				completedStrataQueinnec[curStrata] = true;
				curStrata++;
				continue;
			}

			Index index = strataIndexes->at(nextIndex);
			nextIndexes[curStrata]++;

			double x = GT[0] + index.x * GT[1] + index.y * GT[2];
			double y = GT[3] + index.x * GT[4] + index.y * GT[5];
			OGRPoint newPoint = OGRPoint(x, y);
	
			if (mindist != 0.0 && p_layer->GetFeatureCount() != 0) {
				bool add = true;
				for (const auto &p_feature : *p_layer) {
					OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
					if (newPoint.Distance(p_point) < mindist) {
						add = false;
						break;
					}
				}
			
				if (!add) {
					curStrata++;
					continue;
				}
	
			}

			addPoint(&newPoint, p_layer);
			addedSamples++;
			samplesAddedPerStrata[curStrata]++;

			if (plot) {
				xCoords.push_back(x);
				yCoords.push_back(y);
			}

			curStrata++;
		}
	}
		
	strataIndexVectors = indices.getStrataIndexVectors(samplesAddedPerStrata, strataSampleCounts, rng);	
	for (int64_t i = 0; i < numStrata; i++) {
		//set next indexes to 0 because the vectors are different than the queinnec vectors
		nextIndexes[i] = 0;
	}
	curStrata = 0;
	
	//step 8: generate coordinate points for each sample index.
	while (numCompletedStrata < numStrata && addedSamples < numSamples) {
		if (curStrata == numStrata) {
			curStrata = 0;
		}
		if (completedStrata[curStrata]) {
			curStrata++;
			continue;
		}

		int64_t sampleCount = strataSampleCounts[curStrata];
		int64_t samplesAdded = samplesAddedPerStrata[curStrata];
		if (samplesAdded == sampleCount) {
			numCompletedStrata++;
			completedStrata[curStrata] = true;
			curStrata++;
			continue;
		}

		std::vector<Index> *strataIndexes = strataIndexVectors[curStrata];
		size_t nextIndex = nextIndexes[curStrata];
		if (strataIndexes->size() == nextIndex) {
			numCompletedStrata++;
			completedStrata[curStrata] = true;
			curStrata++;
			continue;
		}

		Index index = strataIndexes->at(nextIndex);
		nextIndexes[curStrata]++;

		double x = GT[0] + index.x * GT[1] + index.y * GT[2];
		double y = GT[3] + index.x * GT[4] + index.y * GT[5];
		OGRPoint newPoint = OGRPoint(x, y);

		if (mindist != 0.0 && p_layer->GetFeatureCount() != 0) {
			bool add = true;
			for (const auto &p_feature : *p_layer) {
				OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
				if (newPoint.Distance(p_point) < mindist) {
					add = false;
					break;
				}
			}
		
			if (!add) {
				curStrata++;
				continue;
			}

		}

		addPoint(&newPoint, p_layer);
		addedSamples++;
		samplesAddedPerStrata[curStrata]++;

		if (plot) {
			xCoords.push_back(x);
			yCoords.push_back(y);
		}

		curStrata++;
	}
	
	//step 9: create GDALVectorWrapper to store dataset of sample points
	GDALVectorWrapper *p_vector = new GDALVectorWrapper(p_samples);

	//step 10: write vector if filename is not "".
	//
	//TODO rather than first making an in-memory dataset then writing to a file afterwards,
	//just make the correct type of dataset from the get go
	if (filename != "") {
		try {
			p_vector->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	size_t actualSampleCount = static_cast<size_t>(p_layer->GetFeatureCount());
	return {{xCoords, yCoords}, p_vector, actualSampleCount};
}

PYBIND11_MODULE(strat, m) {
	m.def("strat_cpp", &strat,
		pybind11::arg("p_raster"),
		pybind11::arg("bandNum"),
		pybind11::arg("numSamples"),
		pybind11::arg("numStrata"),
		pybind11::arg("allocation"),
		pybind11::arg("weights"),
		pybind11::arg("p_mraster").none(true),
		pybind11::arg("mrastBandNum"),
		pybind11::arg("method"),
		pybind11::arg("wrow"),
		pybind11::arg("wcol"),
		pybind11::arg("mindist"),
		pybind11::arg("p_existing").none(true),
		pybind11::arg("force"),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("filename"),
		pybind11::arg("tempFolder"),
		pybind11::arg("threads"));
}
