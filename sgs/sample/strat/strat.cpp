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
		
		this->indexCountPerStrata.resize(numStrata);
		this->indexesPerStrata.resize(numStrata);
		
		this->firstXIndexCountPerStrata.resize(numStrata);
		this->firstXIndexesPerStrata.resize(numStrata);
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
			std::vector<Index>().swap(firstXIndexesPerStrata[i]);
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
			xso::xoshiro_4x64_plus rng = *p_rng;
			randVals[i] = ((rng() >> 11) & multiplier) == multiplier;
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
class FocalWindow {
private:
	int fwHeight, fwWidth;
	std::vector<bool> matrix;
	std::vector<bool> prevVertSame;
	bool prevHoriSame;
	int64_t fwyStart, fwyMidStart, fwyEnd, fwyMidEnd;
	int64_t fwxStart, fwxMidStart, fwxEnd, fwxMidEnd;
	bool blockTopPad, blockBotPad, blockLeftPad, blockRightPad;

public:
	bool scanlines;
	int64_t fwx, fwy;
	int horizontalPad, verticalPad;
	int wrow, wcol;
	int xOff, yOff;
	bool addSelf;
	bool addFw;
	void *p_fwScanline;

	/**
	 *
	 */
	inline void
	init(int width, int xBlockSize, int yBlockSize, int wrow, int wcol, size_t bandPixelSize) {
		this->scanlines = xBlockSize == width;
		this->horizontalPad = wcol / 2;
		this->verticalPad = wrow / 2;
		this->fwHeight = wrow;
		this->fwWidth = xBlockSize - wcol + 1;
		this->matrix.resize(fwWidth * fwHeight, true);
		this->prevVertSame.resize(xBlockSize, true);
		this->prevHoriSame = true;
		this->fwy = -wrow + 1;
		this->p_fwScanline = this->scanlines ? VSIMalloc2(width, bandPixelSize) : nullptr;
		this->wrow = wrow;
		this->wcol = wcol;
	}

	/**
	 *
	 */
	inline void
	readNewBlock(RasterBandMetaData& band, int width, int height, int xBlock, int yBlock, int xBlocks, int yBlocks, int xBlockSize, int yBlockSize, int& xValid, int& yValid) {
		if (!this->scanlines) {
			throw std::runtime_error("SHOULD NOT BE HERE.");
		}
		
		this->blockTopPad = yBlock != 0;
		this->blockBotPad = yBlock != yBlocks - 1;
		this->blockLeftPad = xBlock != 0;
		this->blockRightPad = xBlock != xBlocks - 1;

		this->xOff = xBlock * xBlockSize - blockLeftPad * horizontalPad;
		xValid = std::min(xBlockSize + blockLeftPad * horizontalPad + blockRightPad * horizontalPad, width - xOff);
		this->yOff = yBlock * yBlockSize - blockTopPad * verticalPad;
		yValid = std::min(yBlockSize + blockTopPad * verticalPad + blockBotPad * verticalPad, height - yOff);

		band.p_band->RasterIO(
			GF_Read,
			xOff,
			yOff,
			xValid,
			yValid,
			band.p_buffer,
			xValid,
			yValid,
			band.type,
			0,
			0
		);
	}

	/**
	 *
	 */
	inline void
	updateValsForNewBlock(int xValid) {
		this->fwHeight = this->wrow;
		this->fwWidth = xValid - this->wcol + 1;
		this->matrix.resize(this->fwHeight * fwWidth, true);
		this->prevVertSame.resize(xValid, true);
		this->fwy = -this->wrow + 1;
	}

	/**
	 *
	 */
	inline void
	readInFocalWindowScanline(RasterBandMetaData& band, int width, int yBlock, int yBlockSize, int y) {
		band.p_band->RasterIO(
			GF_Read,
			0,
			yBlock * yBlockSize + y - wrow,
			width,
			1,
			this->p_fwScanline,
			width, 
			1,
			band.type,
			0,
			0
		);
	}

	/**
	 *
	 */
	inline void
	resetOldFocalWindowMatrixSections(void) {
		for (int64_t fwxi = 0; fwxi < this->fwWidth; fwxi++) {
			this->matrix[((fwy + wrow - 1) % wrow) * fwWidth + fwxi] = true;
		}
	}

	/**
	 *
	 */
	inline void
	resetFocalWindowYVals(int height, int yBlock, int yBlockSize, int yValid, int y) {
		this->fwyStart = std::max(this->fwy, static_cast<int64_t>(0));
		this->fwyMidStart = std::max(this->fwy + 1, static_cast<int64_t>(0));
		this->fwyEnd = this->scanlines ? 
			std::min(y + yBlock * yBlockSize, height - this->wrow + 1) :
			std::min(y, yValid - this->wrow + 1);
		this->fwyMidEnd = scanlines ? 
			this->fwyEnd - 1 + static_cast<int64_t>(y + yBlock * yBlockSize > height - this->wrow + 1) :
			this->fwyEnd - 1 + static_cast<int64_t>(y > yValid - this->wrow + 1);
	}

	/**
	 *
	 */
	inline void
	resetfwx(void) {
		this->fwx = -this->wcol + 1;
	}

	/**
	 *
	 */
	inline void
	setAddSelf(int xValid, int yValid, int x, int y) {
		this->addSelf = true; 

		if (!this->scanlines) {
			this->addSelf &= !(this->blockLeftPad && x < horizontalPad) &&
				   	 !(this->blockRightPad && x >= xValid - horizontalPad) &&
				   	 !(this->blockTopPad && y < verticalPad) &&
				   	 !(this->blockBotPad && y >= yValid - verticalPad);

		}
	}

	/**
	 *
	 */
	inline void
	setAddFw(void) {
		this->addFw = this->fwy >= 0 && this->fwx >= 0;
	}

	/**
	 *
	 */
	inline void
	resetFocalWindowXVals(int height, int yBlock, int yBlockSize, int yValid, int x, int y) {
		this->fwxStart = std::max(this->fwx, static_cast<int64_t>(0));
		this->fwxMidStart = std::max(this->fwx + 1, static_cast<int64_t>(0));
		this->fwxEnd = std::min(x, this->fwWidth);
		this->fwxMidEnd = this->scanlines ?
			this->fwxEnd - 1 + static_cast<int64_t>(y + yBlock * yBlockSize > height - this->wrow + 1) :
			this->fwxEnd - 1 + static_cast<int64_t>(y > yValid - this->wrow + 1);
	}

	/**
	 *
	 */
	inline void
	updateFocalWindowMatrix(int x, int y, int yBlock, int yBlockSize, bool nextHoriSame, bool nextVertSame, bool isNan, bool accessible, bool alreadySampled) {
		int64_t fwi;
		//set the portion of the focal window matrix to false which is impacted
		//only when the next horizontal pixel value is different than the 
		//current one.
		if (!nextHoriSame) {
			for (int64_t fwyi = fwyStart; fwyi <= fwyMidEnd; fwyi++) {
				fwi = (fwyi % wrow) * fwWidth + fwxEnd;
				matrix[fwi] = false;
			}
		}

		//set the portion of the focal window matrix to false which is impacted
		//only when the next vertical pixel value is different than the current
		//one.
		if (!nextVertSame) {
			for (int64_t fwxi = fwxStart; fwxi <= fwxMidEnd; fwxi++) {
		 		fwi = (fwyEnd % wrow) * fwWidth + fwxi;
				matrix[fwi] = false;
			}
		}

		//set the portion of the focal window matrix to false which is impacted
		//when either the next vertical or the next horizontal pixel is different
		//than the current one.
		if (!nextHoriSame || !nextVertSame) {
			fwi = (fwyEnd % wrow) * fwWidth + fwxEnd;
			matrix[fwi] = false;
		}

		//set the portion of the focal window matrix to false which must be
		//changed when the next horizontal pixel is different but the previous
		//horizontal pixel is the same as the current one.
		if (!nextHoriSame && prevHoriSame) {
	  		for (int64_t fwxi = fwxMidStart; fwxi <= fwxMidEnd; fwxi++) {
				fwi = (fwyStart % wrow) * fwWidth + fwxi;
				matrix[fwi] = false;
			}
		}

		//set the portion of the focal window matrix to false which must be
		//changed when the next vertical pixel is different but the previous
		//vertical pixel is the same as the current one. 
		if (!nextVertSame && prevVertSame[x]) {
			for (int64_t fwyi = fwyMidStart; fwyi <= fwyMidEnd; fwyi++) {
				fwi = (fwyi % wrow) * fwWidth + fwxStart;
				matrix[fwi] = false;
			}
		}

		//set the portion of the focal window matrix to false which must be
		//changed when either the next vertical or horizontlal pixel is different,
		//but both the previous vertical and previous horizontal pixels were the same
		//as the current one.
		if ((!nextHoriSame || !nextVertSame) && prevHoriSame && prevVertSame[x]) {
			for (int64_t fwyi = fwyMidStart; fwyi <= fwyMidEnd; fwyi++) {
				for (int64_t fwxi = fwxMidStart; fwxi <= fwxMidEnd; fwxi++) {
					fwi = (fwyi % wrow) * fwWidth + fwxi;
					matrix[fwi] = false;
				}
			}
		}

		//if the current pixel is nan, or not accessible, or is already a pre-existing
		//sample, mark it's location in the focal window matrix as false so we know not 
		//to add it in the future without checking for nan or accessibility or existing again.
		if ((isNan || !accessible || alreadySampled) && x - horizontalPad >= 0 && y - verticalPad >= 0) {
			int64_t fwyi = ((scanlines ? y + yBlock * yBlockSize : y) - verticalPad) % wrow;
			int64_t fwxi = x - horizontalPad;
			matrix[fwyi * fwWidth + fwxi] = false;
		}

		this->prevHoriSame = nextHoriSame;
		this->prevVertSame[x] = nextVertSame;
	}

	/**
	 *
	 */
	inline bool
	checkFocalWindowMatrix(void) {
		int64_t fwIndex = (fwy % wrow) * fwWidth + fwx;	
		return matrix[fwIndex];
	}
};

/**
 *
 */
void processBlocksStratRandom(
	RasterBandMetaData& band,
	Access& access,
	Existing& existing,
	RandValController& rand,
	IndexStorageVectors& indices,
	std::vector<int64_t>& existingSampleStrata,
	int width,
	int height) 
{
	int nanInt = static_cast<int>(band.nan);
	int xBlockSize = band.xBlockSize;
       	int yBlockSize = band.yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, band.size);
	int8_t *p_access = nullptr;
	if (access.used) {
		access.band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, access.band.size);
		p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);
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

			//calculate rand vals
			rand.calculateRandValues();

			//iterate through block and update vectors
			for (int y = 0; y < yValid; y++) {
				size_t blockIndex = static_cast<size_t>(y * xBlockSize);
				for (int x = 0; x < xValid; x++) {
					int val = getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex);
					Index index = {x + xBlock * xBlockSize, y + yBlock * yBlockSize};

					bool isNan = val == nanInt;
					bool accessible = !access.used || p_access[blockIndex] == 1;
					bool alreadySampled = existing.used && existing.containsIndex(index.x, index.y);
					blockIndex++;

					//check nan
					if (!isNan) {
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

void processBlocksStratQueinnec(
	RasterBandMetaData &band,
	Access& access,
	Existing& existing,
	RandValController& rand,
	IndexStorageVectors& indices,
	RandValController& queinnecRand,
	IndexStorageVectors& queinnecIndices,
	FocalWindow& focalWindow,
	std::vector<int64_t>& existingSampleStrata,
	int width,
	int height,
	int wrow,
	int wcol) 
{
	int nanInt = static_cast<int>(band.nan);
	int xBlockSize = band.xBlockSize;
	int yBlockSize = band.yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	//allocate buffers
	int8_t *p_access = nullptr;
	if (!focalWindow.scanlines) {
		band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, band.size);
		if (access.used) {
			access.band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, access.band.size);
			p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);
		}	
	}
	else {
		int newXDimSize = xBlockSize + focalWindow.horizontalPad * 2;
		int newYDimSize = yBlockSize + focalWindow.verticalPad * 2;
		band.p_buffer = VSIMalloc3(newXDimSize, newYDimSize, band.size);
		if (access.used) {
			access.band.p_buffer = VSIMalloc3(newXDimSize, newYDimSize, access.band.size);
			p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);
		}
	}

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;

			//read new block
			if (focalWindow.scanlines) {
				band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
				rasterBandIO(band, band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
			}
			else {
				focalWindow.readNewBlock(band, width, height, xBlock, yBlock, xBlocks, yBlocks, xBlockSize, yBlockSize, xValid, yValid);
			}

			//read access block if necessary
			if (access.used) {
				if (focalWindow.scanlines) {
					rasterBandIO(access.band, access.band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
				}
				else {
					access.band.p_band->RasterIO(
						GF_Read,
						focalWindow.xOff,
						focalWindow.yOff,
						xValid,
						yValid,
						access.band.p_buffer,
						xValid,
						yValid,
						access.band.type,
						0,
						0
					);
				}
			}

			//calculate rand values
			rand.calculateRandValues();
			queinnecRand.calculateRandValues();

			if (!focalWindow.scanlines) {
				focalWindow.updateValsForNewBlock(xValid);
			}

			for (int y = 0; y < yValid; y++) {
				if (focalWindow.scanlines && (yBlock * yBlockSize + y - wrow) >= 0) {
					focalWindow.readInFocalWindowScanline(band, width, yBlock, yBlockSize, y);
				}

				size_t blockIndex = focalWindow.scanlines ?
					static_cast<size_t>(y * xBlockSize) :
					static_cast<size_t>(y * yValid);

				focalWindow.resetOldFocalWindowMatrixSections();
				focalWindow.resetFocalWindowYVals(height, yBlock, yBlockSize, yValid, y);
				focalWindow.resetfwx();	

				for (int x = 0; x < xValid; x++) {
					focalWindow.setAddSelf(xValid, yValid, x, y);
					focalWindow.setAddFw();
					focalWindow.resetFocalWindowXVals(height, yBlock, yBlockSize, yValid, x, y);

					int val = getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex);
					Index index = {x + xBlock + xBlockSize, y + yBlock + yBlockSize};

					bool isNan = val == nanInt;
					bool accessible = !access.used || (access.used && p_access[blockIndex] == 1);
					bool alreadySampled = existing.used && existing.containsIndex(index.x, index.y);
					 
					bool nextVertSame = !isNan && ((y == yValid - 1) || 
						val == getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex + (focalWindow.scanlines ? xBlockSize : xValid)));
					bool nextHoriSame = !isNan && ((x == xValid - 1) ||
						val == getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex + 1));

					//UPDATE FOCAL WINDOW MATRIX
					focalWindow.updateFocalWindowMatrix(x, y, yBlock, yBlockSize, nextHoriSame, nextVertSame, isNan, accessible, alreadySampled);
						
					if (focalWindow.addSelf && !isNan) {
						indices.updateStrataCounts(val);

						if (alreadySampled) {
							existingSampleStrata[val]++;
						}

						if (accessible && !alreadySampled) {
							indices.updateFirstXIndexesVector(val, index);	

							if (rand.next()) {
								indices.updateIndexesVector(val, index);
							}
						}
					}

					//ADD NEXT FOCAL WINDOW PIXEL
					if (focalWindow.addFw && focalWindow.checkFocalWindowMatrix()) {
						int val = focalWindow.scanlines ?
								getPixelValueDependingOnType<int>(band.type, focalWindow.p_fwScanline, x + xBlock * xBlockSize - wcol) :
								getPixelValueDependingOnType<int>(band.type, band.p_buffer, (y - wrow) * xValid + (x - wcol));
						index = {x + xBlock * xBlockSize - wcol, y + yBlock * yBlockSize -wrow};
						
						queinnecIndices.updateStrataCounts(val);
						queinnecIndices.updateFirstXIndexesVector(val, index);
						if (queinnecRand.next()) {
							queinnecIndices.updateIndexesVector(val, index);
						}
					}

					focalWindow.fwx++;
				}
				focalWindow.fwy++;
			}
		}
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

	//normal rand used with both methods, queinnec rand used with only queinnec method
	RandValController rand(band.xBlockSize, band.yBlockSize, multiplier, &rng);
	RandValController queinnecRand(band.xBlockSize, band.yBlockSize, queinnecMultiplier, &rng);

	//normal indices used with both methods, queinnec indices used only with queinnec method
	IndexStorageVectors indices(numStrata, 10000);
	IndexStorageVectors queinnecIndices(numStrata, 10000);	

	FocalWindow focalWindow;
	if (method == "Queinnec") {
		focalWindow.init(width, band.xBlockSize, band.yBlockSize, wrow, wcol, band.size);
	}

	if (method == "random") {
		processBlocksStratRandom(
			band,
			access,
			existing,
			rand,
			indices,
			existingSampleStrata,
			width,
			height
		);
	}
	else { //method == queinnec
		processBlocksStratQueinnec(
			band,
			access,
			existing,
			rand,
			indices,
			queinnecRand,
			queinnecIndices,
			focalWindow,
			existingSampleStrata,
			width,
			height,
			wrow,
			wcol
		);
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
