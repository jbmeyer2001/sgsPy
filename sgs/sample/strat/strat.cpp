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
template <typename T>
std::vector<T>
calculateAllocation(
	T numSamples,
	std::string allocation, 
	std::vector<T>& strataCounts, 
	std::vector<double>& weights,
	T numPixels)
{
	std::vector<T> retval;
	T remainder = numSamples;
	T numStrata = strataCounts.size();
	if (allocation == "prop") {
		//allocate the samples per stratum according to stratum size
		T pixelsPerSample = numPixels / numSamples;

		//add 1 if pixelsPerSample was truncated down, to avoid adding too many samples
		pixelsPerSample += static_cast<T>(pixelsPerSample * numSamples < numPixels);

		for (T i = 0; i < numStrata; i++) {
			T count = strataCounts[i] / pixelsPerSample;
			retval.push_back(count);
			remainder -= count;
		}
	}
	else if (allocation == "equal") {
		//determine the count of samples per strata
		T strataSampleCount = numSamples / numStrata;
		
		for (T i = 0; i < numStrata; i++) {
			retval.push_back(strataSampleCount);
			remainder -= strataSampleCount;
		}
	}
	else if (allocation == "manual") {
		//allocate samples accordign to weights.
		for (T i = 0; i < numStrata; i++) {
			T count = static_cast<T>(static_cast<double>(numSamples) * weights[i]);
			retval.push_back(count);
			remainder -= count;
		}
	}
	else { //allocation == "optim"
		//TODO implement
		throw std::runtime_error("'optim' has not been implemented!");
	}

	//redistribute remainder pixels among strata, and check strata sizes
	for (T i = numStrata; i > 0; i--) {
		T extra = remainder / i;
		retval[i - 1] += extra;
		remainder -= extra;

		if (retval[i - 1] > strataCounts[i - 1]) {
			std::cout << "warning: strata " << i - 1 << " does not have enough pixels for the full " << retval[i - 1] << " samples it should recieve. There will be less than " << numSamples << " final samples." << std::endl;
			retval[i - 1] = strataCounts[i - 1];
		}
	}

	return retval;
}

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
template <typename T>
void processBlocksStratRandom(
	RasterBandMetaData& band,
	Access& access,
	Existing& existing,
	IndexStorageVectors& indices,
	std::vector<int64_t>& existingSampleStrata,
	uint64_t multiplier,
	xso::xoshiro_4x64_plus& rng,
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
					blockIndex++;

					//check nan
					if (isNan) {
						continue;
					}

					//udpate strata counts
					indices.updateStrataCounts(val);
					
					//update existing sampled strata
					if (alreadySampled) {
						existingSampleStrata[val]++;
					}

					//add val to stored indices
					if (accessible && !alreadySampled) {
						indices.updateFirstXIndexesVector(val, index);

						if (rand.next()) {
							indices.updateIndexesVector(val, index);
						}					
					}
				}
			}
		}
	}

	VSIFree(band.p_buffer);
	if (access.used) {
		VSIFree(access.band.p_buffer);
	}	
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
void processBlocksStratQueinnec(
	RasterBandMetaData &band,
	Access& access,
	Existing& existing,
	IndexStorageVectors& indices,
	IndexStorageVectors& queinnecIndices,
	FocalWindow& fw,
	std::vector<int64_t>& existingSampleStrata,
	uint64_t multiplier,
	uint64_t queinnecMultiplier,
	xso::xoshiro_4x64_plus& rng,
	int width,
	int height) 
{
	T nanInt = static_cast<T>(band.nan);
	std::cout << "nanInt: " << static_cast<int>(nanInt) << std::endl;

	//adjust blocks to be a large chunk of scanlines
	int xBlockSize = band.xBlockSize;
	int yBlockSize = band.yBlockSize;
	if (xBlockSize != width) {
		xBlockSize = width;
	}
	else {
		yBlockSize = 128;
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

	//get number of blocks
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	//TO REMOVE LATER
	int wrow = fw.wrow;
	int wcol = fw.wcol;

	int fwyi = 0;
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		//read block
		int xOff = 0;
		int yOff = yBlock * yBlockSize;
		int xValid = width;
		int yValid = std::min(yBlockSize, height - yBlock * yBlockSize);

		if (yBlock == 0) { 
			band.p_band->RasterIO(GF_Read, xOff, yOff, xValid, yValid, band.p_buffer, xValid, yValid, band.type, 0, 0);
		}
		else {
			int yOffRead = yOff - fw.vpad * 2;
			int yValidRead = yValid + fw.vpad * 2;
			band.p_band->RasterIO(GF_Read, xOff, yOffRead, xValid, yValidRead, band.p_buffer, xValid, yValidRead, band.type, 0, 0);
		}

		//read access block
		if (access.used) {
			int axOff = 0;
			int ayOff = yBlock * yBlockSize;
			int axValid = width;
			int ayValid = std::min(yBlockSize, height - yBlock * yBlockSize);
			access.band.p_band->RasterIO(GF_Read, axOff, ayOff, axValid, ayValid, access.band.p_buffer, axValid, ayValid, band.type, 0, 0);
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

				//update strata counts
				indices.updateStrataCounts(val);

				//update existing samled strata
				if (alreadySampled) {
					existingSampleStrata[val]++;
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

				//update strata counts
				indices.updateStrataCounts(val);

				//update exisitng sampled strata
				if (alreadySampled) {
					existingSampleStrata[val]++;
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

				//update strata counts
				indices.updateStrataCounts(val);

				//update existing samled strata
				if (alreadySampled) {
					existingSampleStrata[val]++;
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

		std::cout << "completed yBlock = " << yBlock << "/" << yBlocks << std::endl;
	}

	VSIFree(band.p_buffer);
	if (access.used) {
		VSIFree(access.band.p_buffer);
	}	
}

/**
 * This function conducts stratified random sampling on the provided stratified raster.
 *
 *
 * First, the raster is iterated over, the no data pixels, and the pixels which
 * are inaccessible are ignored. The remaining accessable data pixel indexes
 * are placed into vectors corresponding to their strata. After this iteration,
 * there will be a a number of vectors equivalent to the number of strata,
 * and they will contain the indexes of the pixels which are their
 * corresponding strata.
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
	std::string method,
	int wrow,
	int wcol,
	double mindist,
	GDALVectorWrapper *p_existing,
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
	std::vector<int64_t> existingSampleStrata(numStrata, 0);	
	Existing existing(
		p_existing,
		GT,
		width,
		p_layer,
		plot,
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

	if (method == "random") {
		switch (band.type) {
			case GDT_Int8:
				processBlocksStratRandom<int8_t>(band, access, existing,
							 indices, existingSampleStrata,
							 multiplier, rng, width, height);
				break;
			case GDT_Int16:
				processBlocksStratRandom<int16_t>(band, access, existing,
							 indices, existingSampleStrata,
							 multiplier, rng, width, height);
				break;
			default:
				processBlocksStratRandom<int32_t>(band, access, existing,
							 indices, existingSampleStrata,
							 multiplier, rng, width, height);
				break;
		}
			}
	else { //method == queinnec
		switch (band.type) {
			case GDT_Int8:
				processBlocksStratQueinnec<int8_t>(band, access, existing, indices,
							   queinnecIndices, fw, existingSampleStrata,
							   multiplier, queinnecMultiplier, rng, width, height);
				break;
			case GDT_Int16:
				processBlocksStratQueinnec<int16_t>(band, access, existing, indices,
							   queinnecIndices, fw, existingSampleStrata,
							   multiplier, queinnecMultiplier, rng, width, height);
				break;
			default:
				processBlocksStratQueinnec<int32_t>(band, access, existing, indices,
							   queinnecIndices, fw, existingSampleStrata,
							   multiplier, queinnecMultiplier, rng, width, height);
				break;
		}
	}

	std::vector<int64_t> strataCounts = indices.getStrataCounts();
	int64_t numDataPixels = indices.getNumDataPixels();

	std::vector<int64_t> strataSampleCounts = calculateAllocation<int64_t>(
		numSamples,
		allocation,
		strataCounts,
		weights,
		numDataPixels
	);

	std::vector<std::vector<Index> *> strataIndexVectors;
	std::vector<size_t> nextIndexes(numStrata, 0);

	std::vector<bool> completedStrata(numStrata, false);
	std::vector<bool> completedStrataQueinnec(numStrata, false);
	std::vector<int64_t> samplesAddedPerStrata = existingSampleStrata;
 	int64_t numCompletedStrata = 0;
	int64_t numCompletedStrataQueinnec = 0;
	int64_t curStrata = 0;
	int64_t addedSamples = 0;

	for (const int64_t& samples : existingSampleStrata) {
		addedSamples += samples;
	}
	
	if (method == "Queinnec") {
		strataIndexVectors = queinnecIndices.getStrataIndexVectors(samplesAddedPerStrata, strataSampleCounts, rng);	

		size_t curStrata = 0;
		while (numCompletedStrataQueinnec < numStrata && addedSamples < numSamples) {
			if (curStrata == numStrata) {
				curStrata = 0;
			}
			if (completedStrataQueinnec[curStrata]) {
				curStrata++;
				continue;
			}

			size_t sampleCount = strataSampleCounts[curStrata];
			size_t samplesAdded = samplesAddedPerStrata[curStrata];
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
	for (size_t i = 0; i < numStrata; i++) {
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

		size_t sampleCount = strataSampleCounts[curStrata];
		size_t samplesAdded = samplesAddedPerStrata[curStrata];
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
		pybind11::arg("method"),
		pybind11::arg("wrow"),
		pybind11::arg("wcol"),
		pybind11::arg("mindist"),
		pybind11::arg("p_existing").none(true),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("filename"),
		pybind11::arg("tempFolder"));
}
