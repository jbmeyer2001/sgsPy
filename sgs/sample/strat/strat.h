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

#include "utils/access.h"
#include "utils/existing.h"
#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

#include <xoshiro.h>

namespace strat {

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
 * This struct deals with the 'optim' allocation method. The optim allocation
 * method requires that the within-strata variance be calculated based on
 * a seperate raster band. This struct stores band information for the
 * raster band which is used to calculate the variance, the variance per
 * strata, as well as whether the optim method is used.
 */
struct OptimAllocationDataManager {
	RasterBandMetaData band;
	std::vector<Variance> variances;
	bool used = false;

	/**
	 * Constructor, creates an instance of OptimAllocationDataManager given a pointer to the
	 * raster containing the band to use, as well as the integer band number. The allocation
	 * method is also passed, and if the allocation is not 'optim', this struct remains unused.
	 *
	 * @param GDALRasterWrapper *p_raster
	 * @param int bandNum
	 * @param std::string allocation
	 */
	OptimAllocationDataManager(GDALRasterWrapper *p_raster, int bandNum, std::string allocation) {
		if (!p_raster || allocation != "optim") {
			return;
		}

		this->band.p_band = p_raster->getRasterBand(bandNum);
		this->band.type = p_raster->getRasterBandType(bandNum);
		this->band.size = p_raster->getRasterBandTypeSize(bandNum);
		this->band.p_buffer = nullptr; //allocated later
		this->band.nan = this->band.p_band->GetNoDataValue();
		this->band.p_band->GetBlockSize(&this->band.xBlockSize, &this->band.yBlockSize);
		this->used = true;
	}
	
	/**
	 * Deconstructor, this method frees memory which would be allocated if this object
	 * was used.
	 */
	~OptimAllocationDataManager() {
		if (this->used) {
			VSIFree(this->band.p_buffer);
		}
	}

	/**
	 * Initializes some of the data within the struct, setting the size of the vector
	 * which contains Variance information, as well as allocating the memory which
	 * will be used to store the raster band.
	 *
	 * @param int numStrata
	 * @param int xBlockSize
	 * @param int yBlockSize
	 */
	inline void
	init(int numStrata, int xBlockSize, int yBlockSize) {
		this->variances.resize(numStrata);
		this->band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, band.size);
	}

	/**
	 * This function reads a new block of data in from the raster band which this
	 * struct controls.
	 *
	 * @param int xBlockSize
	 * @param int yBlockSize
	 * @param int xBlock
	 * @param int yBlock
	 * @param int xValid
	 * @param int yValid
	 */
	inline void
	readNewBlock(int xBlockSize, int yBlockSize, int xBlock, int yBlock, int xValid, int yValid) {
		rasterBandIO(this->band, this->band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
	}

	/**
	 * This function updates the Variance calculation for a particular pixel. The index of the pixel within
	 * the block is given, to get the pixel value. The strata is also given to indicate which strata's variance
	 * to update. This function should not be called if the pixel is nan.
	 *
	 * @param int index
	 * @param int strata
	 */
	inline void
	update(int index, int strata) {
		double val = getPixelValueDependingOnType<double>(this->band.type, this->band.p_buffer, index);
		variances[strata].update(val);	
	}

	/**
	 * This function gets the allocation percentages per strata using the variances. The optim
	 * allocation method is specified by Gregoire and Valentine https://doi.org/10.1201/9780203498880
	 * Section 5.4.4. The method requires that the count and standard deviation of a particular strata
	 * be multiplied, and the proportion of samples to go to the particular strata strata is this 
	 * product divided by the sum of the product of standard deviation and count for every strata.
	 *
	 * The proportions are calculated using the vector of Variance calculations which have been updated
	 * throughout the iteration if the raster, and returned as a vector of double.
	 *
	 * @returns std::vector<double>
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

/**
 * This is a helper function used for determining the probability multiplier for a given raster.
 * The probability of any given pixel being added is the number of samples divided by the
 * number of total pixels. 
 *
 * Rather than storing the indexes of all possible (accessible, not nan) pixels, which is potentially
 * encredibly memory-inefficient for large rasters, it is much better to only store the indexes
 * of roughly the number of total pixels we need to sample. A random number generator is used 
 * for each pixel which is a candidate for being added as a sample. The sample is retained if 
 * the random number generator creates a series of  all 1's for the first n bits. This value n determines 
 * the probability a sample is added. For example, if n were three then 1/8 or 1/(2^n) proportoin of pixels 
 * would be sampled.
 *
 * It can be seen that setting up an n value which is close to the probability samples/pixels but
 * an over estimation would result in an adequte number of indexes stored WITHOUT storing a
 * rediculous number of values.
 *
 * The way this number n is enforced, is by determining a multiplier that takes the form of the first
 * n bits are 1 and the remaining are 0. For example:
 * 1 	-> 00000001 	-> 50%
 * 3 	-> 00000011 	-> 25%
 * 7 	-> 00000111 	-> 12.5%
 * 63 	-> 00111111	-> 1.56%
 *
 * The AND of this multiplier is taken with the rng value, and the result of that and is compared against
 * the multiplier. The 0's from the and remove the unimportant bits, and the 1's enforce the first n
 * values at the beginning.
 *
 * The multiplier is determined by determining the numerator and denominator of this probability (samples/pixels),
 * with extra multipliers for an unknonwn amount of nan values, and multiplying by extra if the mindist parameter
 * is passed as it may cause samples to be thrown out. Further, if an access vector is given and all samples 
 * must fall within the accessible area, the probability is increased by the ratio of the total area in the raster
 * to the accessible area. The probability would then simply be numerator/denominator, but we want a multiplier
 * with a specific number of bits not a small floating point value. The log base 2 is used to transform this division
 * int a subtraction problem, resulting in the number of bits. The value 1 is then left shifted by the number of bits,
 * and subtracted by 1 to give the multiplier.
 *
 * @param GDALRasterWrapper *p_raster
 * @param int numSamples
 * @param bool useMindist
 * @param double accessibleArea
 * @param bool queinnec
 *
 * @returns uint64_t
 */
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
 * This struct is responsible for storing the indices saved while iterating through the raster. 
 *
 * The strataCounts vector strores the total count of each strata within the stratified raster. 
 * The indexCountPerStrata vector stores the total number of indices saved for each strata.
 * The indexesPerStrata vector stores the indices saved of each strata.
 * 
 * Since the pixels are probabilistically added, and this probability is calculated before
 * knowing the total proportion of strata within the strat raster, there is a potential problem
 * if one strata does not appear as much as the others, as this could mean less than the
 * total number of required pixels for that strata are actually saved. This is where the
 * 'first x' vectors come in. They always store the first x number of pixels in a particular
 * strata, meaning that if there are less than x (10,000) pixels of a particular strata,
 * all pixels within that strata will be saved and there is not a worry of missing out on
 * a particular strata due to having a lower proportion. 
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
	 * Constructor. Set the numStrata and x variables.
	 * Also, resize the vectors so that their size is numStrata.
	 *
	 * @param int64_t numStrata
	 * @param int64_t x
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
	 * update the total count by 1 of a particular strata.
	 *
	 * @param int strata
	 */
	inline void
	updateStrataCounts(int strata) {
		this->strataCounts[strata]++;
	}

	/**
	 * update the indexes vector with a new index of a particular
	 * strata.
	 *
	 * @param int strata
	 * @param Index& index
	 */
	inline void
	updateIndexesVector(int strata, Index& index) {
		this->indexesPerStrata[strata].push_back(index);
		this->indexCountPerStrata[strata]++;
	}

	/**
	 * update the first x indexes vector with a new index
	 * of a particular strata.
	 *
	 * @param int strata
	 * @param Index& index
	 */
	inline void
	updateFirstXIndexesVector(int strata, Index& index) {
		int64_t i = firstXIndexCountPerStrata[strata];
		if (i < this->x) {
			firstXIndexesPerStrata[strata][i] = index;
			firstXIndexCountPerStrata[strata]++;	
		}
		else if (i == this->x) {
			std::vector<Index>().swap(firstXIndexesPerStrata[strata]);
			firstXIndexCountPerStrata[strata]++;
		}
	}

	/**
	 * get references to shuffled strata index vectors. This method is called after the 
	 * whole raster has been iterated through, and potential sample pixels are saved.
	 *
	 * If there are not enough values in the normal saved strata, but the first X samples
	 * have been saved (and there hasn't been too many pixels), then the first X samples
	 * are used. Otherwise, the normal saved samples are used.
	 *
	 * A reference is saved for a shuffled version of the vector chosen. The reason a shuffled
	 * version is returned is because the raster is iterated through in an ordered manor (by blocks).
	 * If it wasn't shuffled, the samples selected would be the ones closes to the beginning.
	 *
	 * @param std::vector<int64_t> existing
	 * @param std::vector<int64_t> strataSampleCounts
	 * @param xso::xoshiro_4x64_plus& rng
	 *
	 * @returns std::vector<std::vector<Index> *>
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
	 * get the pixel counts per strata.
	 *
	 * @returns std::vector<int64_t>
	 */
	inline std::vector<int64_t>
	getStrataCounts(void) {
		return this->strataCounts;
	}

	/**
	 * get the total number of data (not nan) pixels.
	 *
	 * @returns int64_t
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
 * This struct controls the calculation and usage of random values during the
 * iteration through the raster. A random number must be generated for
 * each pixel to see if it will be saved for potential sampling.
 *
 * The xoshiro random number generator is used because it is efficient and 
 * statistically sound. The specific generator used (xso::xoshrio_4x64_plus) is used
 * because it is very fast. However, it's lowest 11 bits have low linear complexity (Blackman & Vigna).
 * 
 * We have no need for these lower 11 bits, instead using only the upper 53 bits of the uint64_t value.
 * Proof of this is that, supposing we require the use of all 53 bits, this means a probability of 
 * 1/(2^(56)), or roughly 1 sample per 10^16 pixels. If there were 10^16 pixels to process than a minimum of
 * multiple years would likely pass before execution finished.
 *
 * Rather than calling the generator on every iteration, the generator is repeatedly called at the 
 * beginning of a block for the remaining required pixels, and the true/false values for whether
 * to save a pixel or not are stored in a vector of type boolean.
 */
class RandValController {
private:
	std::vector<bool> randVals;
	size_t randValIndex = 0;
	uint64_t multiplier = 0;
	xso::xoshiro_4x64_plus *p_rng = nullptr;

public:
	/**
	 * Constructor, sets the size of the boolean vector, and assigns the randValIndex, multiplier, and p_rng
	 * member variables.
	 *
	 * @param int xBlockSize
	 * @param int yBlockSize
	 * @param uint64_t multiplier
	 * @param xso::xoshiro_4x64_plus *p_rng
	 */
	RandValController(int xBlockSize, int yBlockSize, uint64_t multiplier, xso::xoshiro_4x64_plus *p_rng) {
		this->randVals.resize(xBlockSize * yBlockSize);
		this->randValIndex = static_cast<size_t>(xBlockSize * yBlockSize);
		this->multiplier = multiplier;
		this->p_rng = p_rng;
	}

	/**
	 * Calculates the true/false values from rand values, a number of times equal to the
	 * number of used random values from the previous block. The return value of the
	 * random number generator is bit shifted by 11 to ignore the lower 11 bits, which
	 * have low linear complexity.
	 *
	 * Next, the bit shifted random value is masked with the multiplier, and if the random
	 * value contains a 1 in every bit which the multiplier does, true is added to the rand
	 * val vector.
	 *
	 * This function is called before iterating through a new block.
	 */
	inline void 
	calculateRandValues(void) {
		for (size_t i = 0; i < randValIndex; i++) {
			randVals[i] = (((*p_rng)() >> 11) & multiplier) == multiplier;
		}
		randValIndex = 0;
	}

	/**
	 * get the next boolean value from the storage vector, and iterate the index to this
	 * vector.
	 */
	inline bool 
	next(void) {
		bool retval = randVals[randValIndex];
		randValIndex++;
		return retval;
	}
};

/**
 * This function processes the strat raster in blocks using the 'random' method. In the random
 * method, every pixel in a particular strata has the same priority of being added as any
 * other pixel in that strata.
 *
 * First, memory is allocated and structures are initialized for random value calculation and
 * optim allocation.
 *
 * Next, iterate through the blocks within the raster. For each block:
 *  - read the strat raster (and potentially access & optim rasters) into memory
 *  - calculate rand values for the new block
 *  - iterate through the pixels in the block
 *
 * For each pixel:
 *  - ignore the pixel if it is a nan value
 *  - update the optim variances if optim allocation used
 *  - update total sample counts
 *  - add index to existing vector if the pixel is part of an existing sample network
 *  - If the pixel is both accessible and not already sampled, update the index storage vectors
 *
 * Once all blocks have been processed, free any allocated memory and calculate the sample
 * allocation per strata, and return this allocation.
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
	int height) 
{
	T nanInt = static_cast<T>(band.nan);
	int xBlockSize = band.xBlockSize;
       	int yBlockSize = band.yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, band.size);
	T *p_buffer = reinterpret_cast<T *>(band.p_buffer);
	int8_t *p_access = nullptr;
	if (access.used) {
		access.band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, access.band.size);
		p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);
	}

	RandValController rand(band.xBlockSize, band.yBlockSize, multiplier, &rng);

	if (optim.used) {
		optim.init(numStrata, xBlockSize, yBlockSize);
	}

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
	
			//read block
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			rasterBandIO(band, band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true, false);

			//read access block
			if (access.used) {
				rasterBandIO(access.band, access.band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
			}

			//read mraster block for optim
			if (optim.used) {
				optim.readNewBlock(xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid);
			}

			//calculate rand vals
			rand.calculateRandValues();

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
					if (optim.used) {
						optim.update(blockIndex, val);
					}

					//udpate strata counts
					indices.updateStrataCounts(val);
					
					//update existing sampled strata
					if (alreadySampled) {
						existingSamples[val].push_back(existing.getPoint(index.x, index.y));
					}

					//add val to stored indices
					if (accessible && !alreadySampled) {
						indices.updateFirstXIndexesVector(val, index);

						if (rand.next()) {
							indices.updateIndexesVector(val, index);
						}					
					}

					//increment block index
					blockIndex++;
				}
			}
		}
	}

	VSIFree(band.p_buffer);
	if (access.used) {
		VSIFree(access.band.p_buffer);
	}	

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
 * The focal window is a method used by Queinnec sampling to ensure that
 * pixels which are surrounded by other pixels of the same strata are prioritized
 * over pixels which aren't.
 *
 * This struct is used to control the storage of values to keep track of which
 * pixels are eligible as 'queinnec' values i.e. surrounded by pixels of the
 * same strata.
 *
 * the m vector is a vector of bool which keeps track whether a particular pixels
 * is eligible considering it's horizontal pixels. The values are initialized to false,
 * and when iterating through the raster, the m vector is set if it is horizontally 
 * surrounded by the same pixels. These saved values will then be used later to determine
 * if a whole pixel is eligible. The logic is used that if all vertical values within the
 * horizontal range of the focal window are the same, AND this stored horizontal value
 * is true for all of the pixels, then the pixel in the middle must be surrounded by
 * whole pixels of the same type.
 *
 * the valid vector is a vector of bool which keeps track whether a particular pixel is
 * available for sampling. Specifically, it isn't nan, is accessible, and isn't an 
 * existing sample. The reason why this data needs to be saved, rather than used for
 * the current pixel, is that in order to check the vertical focal window we must have
 * iterated past a particular pixel before we know whether it is fully eligible. 
 * Rather than checking if it is available for sampling twice, the second time
 * for the focal window's purpose, this value is saved in the valid vector.
 *
 * Both the m and the valid vectors have the same width of the raster, but significantly
 * less height. They only need to have the height of the focal window, and when an
 * old row becomes unused it is reset to false and used for another row in the raster.
 */
struct FocalWindow {
	int wrow, wcol;
	int width;
	int vpad, hpad;
	std::vector<bool> m;
	std::vector<bool> valid;

	/**
	 * Constructor, uses wrow, wcol, and width values to set corrosponding
	 * values. Also calculate vertical pad (vpad), horizontal pad (hpad),
	 * and resize the m and valid vectors.
	 *
	 * @param int wrow
	 * @param int wcol
	 * @param int width
	 */
	FocalWindow(int wrow, int wcol, int width) {
		this->wrow = wrow;
		this->vpad = wrow / 2;
		this->wcol = wcol;
		this->hpad = wcol / 2;
		this->width = width;
		this->m.resize(wrow * width, false);
		this->valid.resize(wrow * width, false);
	}

	/**
	 * this function is used to reset the values to fales for a no longer
	 * used row, so that it can be used for the next row.
	 *
	 * This function is called at the beginning of each scanline/row.
	 *
	 * @param int row
	 */
	inline void
	reset(int row) {
		int start = (row % wrow) * this->width;
		int end = start + this->width;
		for (int i = start; i < end; i++) {
			m[i] = false;
			valid[i] = false;
		}
	}

	/**
	 * This function checks whether a particular index is eligible
	 * to be added as a queinnec pixel, i.e. is surrounded by pixels
	 * of the same strata. 
	 *
	 * Both the m and the valid vectors are checked in this function,
	 * although both the m and valid vectors are set outside of the
	 * functions of this struct. They are set directly during
	 * the iteration through the raster.
	 *
	 * This function only checks whether the pixel is valid (not nan,
	 * accessible, not part of an existing sampling network), and
	 * whether all of the pixels vertically within the focal window 
	 * are horizontally surrounded by the same values. If this
	 * function returns true, those vertical pixels will then be
	 * checked to ensure they are equivalent.
	 *
	 * @param int x
	 * @param int y
	 *
	 * @returns bool
	 */
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
 * This function processes the strat raster in blocks using the 'Queinnec' method. In the Queinnec
 * method, pixels which are surrounded by pixels of the same strata are prioritized for sampling 
 * over pixels which aren't.
 *
 * First, the block size is adjusted to be scanlines with a height of either the original y block size,
 * or 128. This is done because the raster needs to be read as scanlines for the FocalWindow struct
 * to work well and not get any more complicated than it already is. However, we want the raster IO
 * to still be as efficient as possible, so chunks of the raster are still read on block boundaries,
 * just containing a lot more than just 1 block.  
 *
 * Next, memory is allocated and structs are created for random value calculation and optim allocation.
 * More memory is allocated than just 1 of the chunks (xBlockSize * yBlockSize), this is because
 * usign the focal window struct method, we may have to read in some of the final few pixels of the
 * previous chunk, to the start of the new chunk.
 *
 * Next, iterate thorugh the blocks within the raster. For each block:
 *  - read block from the strat raster
 *  - calculate rand values (both normal and queinnec) for the new block
 *  - iterate through the pixels in the block
 *
 * The newBlockStart value is used due to the aforementioned padding which may be placed on the top
 * of the new block. This newBlockStart value is the index of the block which is the start of the
 * new block. This value will be 0 in the case of the first block, but different otherwise.
 *
 * Rather than using a typical nested for loop, one for the vertical and one for the horizontal,
 * we use one for the vertical and three for the horizontal. This is because there are areas on the
 * left and right of the raster which will never be eligible as queinnec pixels because their focal window
 * includes pixels which go off the edge of the raster. And, critically, trying to calculate whether they
 * are horizontally eligible as a queinnec pixel would result in checking a pixel which either doesn't exist
 * or is on the opposite side of the raster.
 *
 * For each pixel horizontally not eligible to be a queinnec pixel:
 *  - ignore the pixel if it is a nan value
 *  - update the optim variances if optim allocation is used
 *  - update total sample counts
 *  - add index to existing vector if the pixel is part of an existing sample network
 *  - if the pixel is both accessible and not already sampled, update the index storage vectors.
 *
 * For pixel which is horizontally eligible to be a queinnec pixel:
 *  - do all of the same as those non-eligible pixels in addition to...
 *  - set focal window valid vector for the current pixel
 *  - set focal window matrix vector for the current pixel by checking horizontally adjacent pixels
 *  - check the focal window matrix for the pixel which will have just had all of it's vertical pixels horizontally checked.
 *    Calling this check on a pixel which would have a negative y will always result in a false due to the focal window
 *    matrix being automatically set to false. If this check succeeds -- check vertical pixels to see if they are the same 
 *    and if so update the queinnec index storage vectors.
 *
 * Once all blocks have been processed, free any allocated memory and calculate the sample allocation per strata,
 * and return this allocation.
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
	FocalWindow& fw,
	std::vector<std::vector<OGRPoint>>& existingSamples,
	uint64_t multiplier,
	uint64_t queinnecMultiplier,
	xso::xoshiro_4x64_plus& rng,
	std::string allocation,
	OptimAllocationDataManager& optim,
	std::vector<double> weights,
	int width,
	int height) 
{
	T nanInt = static_cast<T>(band.nan);

	//adjust blocks to be a large chunk of scanlines
	int xBlockSize = band.xBlockSize;
	int yBlockSize = band.yBlockSize;
	if (xBlockSize != width) {
		xBlockSize = width;
	}
	else {
		yBlockSize = std::min(128, height);
	}

	//allocate required memory
	band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize + fw.vpad * 2, band.size);
	T *p_buffer = reinterpret_cast<T *>(band.p_buffer);
	int8_t *p_access = nullptr;
	if (access.used) {
		access.band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, access.band.size);
		p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);
	}

	RandValController rand(xBlockSize, yBlockSize, multiplier, &rng);
	RandValController queinnecRand(xBlockSize, yBlockSize, queinnecMultiplier, &rng);

	if (optim.used) {
		optim.init(numStrata, xBlockSize, yBlockSize);
	}

	//get number of blocks
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	int fwyi = 0;
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		//read block
		int xOff = 0;
		int yOff = yBlock * yBlockSize;
		int xValid = width;
		int yValid = std::min(yBlockSize, height - yBlock * yBlockSize);

		//read block
		if (yBlock == 0) { 
			band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, band.p_buffer, xValid, yValid, band.type, 0, 0);
		}
		else {
			int stratYOff = yOff - fw.vpad * 2;
			int stratYValid = yValid + fw.vpad * 2;
			band.p_band->RasterIO(GF_Read, xOff, stratYOff, xValid, stratYValid, band.p_buffer, xValid, stratYValid, band.type, 0, 0);
		}

		//read access block
		if (access.used) {
			access.band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, access.band.p_buffer, xValid, yValid, GDT_Int8, 0, 0);
		}

		//read mraster block for optim
		if (optim.used) {
			optim.band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, optim.band.p_buffer, xValid, yValid, optim.band.type, 0, 0);
		}


		//calculate rand vals
		rand.calculateRandValues();
		queinnecRand.calculateRandValues();

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
				if (optim.used) {
					optim.update(y * width + x, val);
				}

				//update strata counts
				indices.updateStrataCounts(val);

				//update existing samled strata
				if (alreadySampled) {
					existingSamples[val].push_back(existing.getPoint(index.x, index.y));
				}

				//add val to stored indices
				if (accessible && !alreadySampled) {
					indices.updateFirstXIndexesVector(val, index);

					if (rand.next()) {
						indices.updateIndexesVector(val, index);
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
				if (optim.used) {
					optim.update(y * width + x, val);
				}

				//update strata counts
				indices.updateStrataCounts(val);

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
					indices.updateFirstXIndexesVector(val, index);

					if (rand.next()) {
						indices.updateIndexesVector(val, index);
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
				if (optim.used) {
					optim.update(y * width + x, val);
				}

				//update strata counts
				indices.updateStrataCounts(val);

				//update existing sampled strata
				if (alreadySampled) {
					existingSamples[val].push_back(existing.getPoint(index.x, index.y));
				}

				//add val to stored indices
				if (accessible && !alreadySampled) {
					indices.updateFirstXIndexesVector(val, index);

					if (rand.next()) {
						indices.updateIndexesVector(val, index);
					}
				}
			}

			fwyi += width;
			if (fwyi == width * fw.wrow) {
				fwyi = 0;
			}
		}
	}

	VSIFree(band.p_buffer);
	if (access.used) {
		VSIFree(access.band.p_buffer);
	}

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
 * First, metadata is acquired in the strat raster band, which is to be read to determine
 * sample strata and to ensure samples don't occur on nan values.
 *
 * Next, the output vector dataset is creates as an in-memory dataset. If the user
 * specifies a filename, this in-memory dataset will be written to disk in a different
 * format after all points have been added.
 *
 * An Access struct is created, which creates a raster dataset containing
 * a rasterized version of access buffers. This raster will be 1 over
 * accessible areas. In the case where there is no access vector given,
 * the structs 'used' member will be false and no processing or rasterization
 * will be done.
 *
 * An Existing struct is created, which retains information on already existing
 * sample points passed in the form of a vector dataset. The points are iterated
 * through and added to the output dataset. The points are also added to a set,
 * and during iteration the indexes of every pixel will be checked against this set
 * to ensure there are no duplicate pixels. In the case whre there is no existing
 * vector given, the structs 'used' member will be false and no processing
 * will be done.
 *
 * Next, a rng() function is created usign the xoshiro library, the specific
 * randm number generator is the xoshrio256++
 * https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf	
 *
 * The impetus behind usign the rng() function to determine which pixels
 * should be added DURING iteration, rather than afterwards, is it removes the
 * necessity of storing every available pixel, which quickly becomes extrordinarily
 * inefficient for large rasters. Rather, for pixels which are accessible, not nan,
 * and not already existing, there is a pre-determined percentage chance to be stored
 * which uses this random number generator. An over-estimation for the percentage
 * chance is made, because it is better to have too many than not enough possible options
 * to sample from. This over-estimation might result in the storage of 2x-3x extra pixels
 * rather than the many orders of magnitude extra storage of adding all pixels. The
 * calculation for this percentage is done and explained in detail in the
 * getProbabilityMultiplier() function.
 *
 * Next, the raster is processed in blocks either using the 'random' or 'Queinnec' 
 * methods, and the return of those functions contains the allocation of samples
 * per strata.
 *
 * Strata are iterated through, with samples being added according to their total allocation.
 * First, existing pixels are added, all of which are added in the case where the force
 * parameter is true. Next, queinnec pixels are added if the queinnec method is used. Finally,
 * remaining random pixels are added.
 *
 * @param GDALRasterWrapper *p_raster
 * @param int bandNum
 * @param int64_t numSamples
 * @param int64_t numStrata
 * @param std::string allocation
 * @param std::vector<double> weights
 * @param GDALRasterWrapper *p_mraster
 * @param int mrasterBandNum
 * @param std::string method
 * @param int wrow
 * @param int wcol
 * @param double mindist
 * @param GDALVectorWrapper *p_existing
 * @param bool force
 * @param GDALVectorWrapper *p_access
 * @param std::string layerName
 * @param double buffInner
 * @param double buffOuter
 * @param bool plot
 * @param std::string filename
 * @param std::string tempFolder
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
	std::string tempFolder)
{
	GDALAllRegister();

	bool useMindist = mindist != 0;
	int width = p_raster->getWidth();
	int height = p_raster->getHeight();
	double *GT = p_raster->getGeotransform();

	std::mutex bandMutex;
	std::mutex rngMutex;
	std::mutex accessMutex;

	//step 1: get raster band
	RasterBandMetaData band;

	GDALRasterBand *p_band = p_raster->getRasterBand(bandNum);
	band.p_band = p_band;
	band.type = p_raster->getRasterBandType(bandNum);
	band.size = p_raster->getRasterBandTypeSize(bandNum);
	band.p_buffer = nullptr;
	band.nan = p_band->GetNoDataValue();
	band.p_mutex = &bandMutex;
	p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);

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

	FocalWindow fw(wrow, wcol, width);

	OptimAllocationDataManager optim(p_mraster, mrastBandNum, allocation);

	std::vector<int64_t> strataSampleCounts; 
	if (method == "random") {
		switch (band.type) {
			case GDT_Int8:
				strataSampleCounts = processBlocksStratRandom<int8_t>(numSamples, numStrata, band, access, 
										      existing, indices, existingSamples,
							 			      multiplier, rng, allocation, optim,
										      weights, width, height);
				break;
			case GDT_Int16:
				strataSampleCounts = processBlocksStratRandom<int16_t>(numSamples, numStrata, band, access, 
										      existing, indices, existingSamples,
							 			      multiplier, rng, allocation, optim,
										      weights, width, height);
				break;
			default:
				strataSampleCounts = processBlocksStratRandom<int32_t>(numSamples, numStrata, band, access, 
										      existing, indices, existingSamples,
							 			      multiplier, rng, allocation, optim,
										      weights, width, height);
				break;
		}
	}
	else { //method == queinnec
		switch (band.type) {
			case GDT_Int8:
				strataSampleCounts = processBlocksStratQueinnec<int8_t>(numSamples, numStrata, band, access, existing, 
											indices, queinnecIndices, fw, existingSamples,
							   				multiplier, queinnecMultiplier, rng, allocation, 
											optim, weights, width, height);
				break;
			case GDT_Int16:
				strataSampleCounts = processBlocksStratQueinnec<int16_t>(numSamples, numStrata, band, access, existing, 
											indices, queinnecIndices, fw, existingSamples,
							   				multiplier, queinnecMultiplier, rng, allocation, 
											optim, weights, width, height);
				break;
			default:
				strataSampleCounts = processBlocksStratQueinnec<int32_t>(numSamples, numStrata, band, access, existing, 
											indices, queinnecIndices, fw, existingSamples,
							   				multiplier, queinnecMultiplier, rng, allocation, 
											optim, weights, width, height);
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

} //namespace strat
