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
 * Next, the calculate_allocation function is used to determine the the total
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
	size_t numSamples,
	size_t numStrata,
	std::string allocation,
	std::vector<double> weights,
	std::string method,
	int wrow,
	int wcol,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string filename,
	bool largeRaster,
	std::string tempFolder)
{
	GDALAllRegister();

	bool useMindist = mindist != 0;
	int width = p_raster->getWidth();
	int height = p_raster->getHeight();

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

	Access access;
	updateAccess(
		p_access, 
		p_raster, 
		layerName, 
		buffInner, 
		buffOuter, 
		largeRaster, 
		tempFolder, 
		band.xBlockSize,
		band.yBlockSize,
		access
	);

	int xBlockSize = band.xBlockSize;
	int yBlockSize = band.yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	band.p_buffer = VSIMalloc3(xBlockSize, yBlockSize, band.size);
	std::vector<bool> randVals(xBlockSize * yBlockSize);
	int randValIndex = xBlockSize * yBlockSize;
	int nanInt = static_cast<int>(band.nan);
	
	uint8_t *p_accessMask = nullptr;
	if (access.used) {
		p_accessMask = static_cast<uint8_t *>(largeRaster ?
			VSIMalloc2(xBlockSize, yBlockSize) :
			access.band.p_buffer);
	}

	std::vector<uint64_t> strataCounts(numStrata, 0);
	std::vector<std::vector<Index>> indexesPerStrata(numStrata);
	std::vector<size_t> indexCountPerStrata(numStrata, 0);

	size_t firstX = 10000;
	std::vector<std::vector<Index>> firstXIndexesPerStrata(numStrata);
	std::vector<size_t> firstXIndexCountPerStrata(numStrata, 0);
	
	//fast random number generator using xoshiro256+
	//https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf
	xso::xoshiro_4x64_plus rng; 

	//the multiplier which will be multiplied by the 53 most significant bits of the output of the
	//random number generator to see whether a pixel should be added or not. The multiplier is
	//a uint64_t number where the least significant n bits are 1 and the remaining are 0. The pixel
	//is added when the least significant n bits (of the bit shifted 53 bits) within the rng match
	//those of the multiplier. The probability a pixel is add is then (1/2^n). Using this method,
	//knowing the amount of pixels, estimating the number of nan pixels, taking into account mindist
	//and accessible area, we can estimate a percentage chance for each pixel and set up a multiplier
	//to make that percentage happen. Doing this enables retaining only a small portion of pixel data
	//and reducing memory footprint significantly, otherwise the index of every pixel in each strata
	//would have to be stored, which would not be feasible for large rasters.
	//
	//NOTE that this method is imprecise, and has an achilles heel where if there are not very many
	//pixels of a particular strata, they may not show up at all in the final samples. For this reason,
	//the first 10000 pixels of every strata are added no matter what. If there aren't enough pixels
	//of a particular strata from the normal probability method, but there are more than 10000 pixels
	//of that strata in the raster, we may have some problems because taking from the 10000 would
	//no longer be random. Perhaps there is a way to dynamically adjust the size of the first 'thousand'
	//pixels as the raster is iterating to take this into account, without retaining too many pixels
	//to be feasible for large images.
	uint64_t multiplier = getProbabilityMultiplier(p_raster, numSamples, useMindist, access.area, false);

	//queinnec method specific vectors
	bool queinnec = method == "Queinnec";
	std::vector<uint64_t> queinnecStrataCounts(numStrata * queinnec, 0);
	std::vector<std::vector<Index>> queinnecIndexesPerStrata(numStrata * queinnec);
	std::vector<size_t> queinnecIndexCountPerStrata(numStrata * queinnec, 0);

	size_t queinnecFirstX = 10000 * queinnec;
	std::vector<std::vector<Index>> queinnecFirstXIndexesPerStrata(numStrata * queinnec);
	std::vector<size_t> queinnecFirstXIndexCountPerStrata(numStrata * queinnec, 0);

	uint64_t multiplierQ = 0;
	if (queinnec) {
		multiplierQ = getProbabilityMultiplier(p_raster, numSamples, useMindist, access.area, true);
	}

	bool scanlines = xBlockSize == width; 
	int horizontalPad = wcol / 2;
	int verticalPad = wrow / 2;
	int fwHeight = wrow;
	int fwWidth = xBlockSize - wcol + 1;
	std::vector<bool> focalWindowMatrix(fwWidth * fwHeight, true);
	std::vector<bool> prevVertSame(xBlockSize, true);
	bool prevHoriSame = true;
	bool nextVertSame;
	bool nextHoriSame;
	int64_t fwy = -wrow + 1;
	int64_t fwx;
	int64_t fwyStart, fwyMidStart, fwyEnd, fwyMidEnd;
	int64_t fwxStart, fwxMidStart, fwxEnd, fwxMidEnd;
	int xOff, yOff;
	bool addSelf;
	bool addFw;
	bool blockTopPad, blockBotPad, blockLeftPad, blockRightPad;
	if (queinnec && !scanlines) {
		band.p_buffer = VSIRealloc(band.p_buffer, (xBlockSize + horizontalPad * 2) * (yBlockSize + verticalPad * 2) * band.size);
	}

	std::vector<bool> queinnecRandVals(xBlockSize * yBlockSize * queinnec);
	int queinnecRandValIndex = xBlockSize * yBlockSize * queinnec;

	void *p_fwScanline = scanlines ? VSIMalloc2(width, band.size) : nullptr;

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {	
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;

			//READ BLOCK
			if (!queinnec || scanlines) {
				//bandMutex.lock();
				band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
				rasterBandIO(band, band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true); 
				//bandMutex.unlock();									//read = true
			}
			else { //queinnec processing by block
				blockTopPad = yBlock != 0;
				blockBotPad = yBlock != yBlocks - 1;
				blockLeftPad = xBlock != 0;
				blockRightPad = xBlock != xBlocks - 1;

				xOff = xBlock * xBlockSize - blockLeftPad * horizontalPad;
				xValid = std::min(xBlockSize + blockLeftPad * horizontalPad + blockRightPad * horizontalPad, width - xOff);
				yOff = yBlock * yBlockSize - blockTopPad * verticalPad;
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

			//READ ACCESS BLOCK IF NECESSARY
			if (access.used) {
				if (!queinnec || scanlines) {
					//accessMutex.lock();
					rasterBandIO(access.band, p_accessMask, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true);
					//accessMutex.lock();									      //read = true
				}
				else {
					access.band.p_band->RasterIO(
						GF_Read,
						xOff,
						yOff,
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
			
			//CALCULATE RAND VALUES
			//rngMutex.lock();
			for (int i = 0; i < randValIndex; i++) {
				randVals[i] = ((rng() >> 11) & multiplier) == multiplier;
			}
			randValIndex = 0;
			//rngMutex.unlock();

			if (queinnec) {
				for (int i = 0; i < queinnecRandValIndex; i++) {
					queinnecRandVals[i] = ((rng() >> 11) & multiplierQ) == multiplierQ;
				}
				queinnecRandValIndex = 0;
			}

			//ITERATE THROUGH BLOCK AND UPDATE VECTORS
			if (!queinnec) {
				for (int y = 0; y < yValid; y++) {
					size_t blockIndex = static_cast<size_t>(y * xBlockSize);
					for (int x = 0; x < xValid; x++) {
						//GET VAL
						int val = getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex);
						
						//CHECK NAN
						bool isNan = val == nanInt;
						if (isNan) {
							blockIndex++;
							continue;
						}
	
						//CHECK ACCESS
						if (access.used && p_accessMask[blockIndex] == 0) {
							strataCounts[val]++;
							blockIndex++;
							continue;
						}
	
						//CREATE INDEX STRUCTURE
						Index index = {x + xBlock * xBlockSize, y + yBlock * yBlockSize};
						strataCounts[val]++;
	
						//UPDATE FIRST X VALS AUTOMATICALLY
						if (firstXIndexCountPerStrata[val] < firstX) {
							int i = firstXIndexCountPerStrata[val];
							firstXIndexCountPerStrata[val]++;
							firstXIndexesPerStrata[val][i] = index;
						}
	
						//CHECK RNG
						if (!randVals[randValIndex]) {
							randValIndex++;
							blockIndex++;
							continue;
						}
						randValIndex++;

						//UPDATE VECTORS
						strataCounts[val]++;
						indexesPerStrata[val].push_back(index);
						indexCountPerStrata[val]++;

						//INCREMENT WITHIN-BLOCK INDEX
						blockIndex++;
					}
				}
			}
			else {
				if (!scanlines) {
					fwHeight = wrow;
					fwWidth = xValid - wcol + 1;
					focalWindowMatrix.resize(fwHeight * fwWidth, true);
					prevVertSame.resize(xValid, true);
					fwy = -wrow + 1;	
				}
					
				for (int y = 0; y < yValid; y++) {
					//IF SCANLINES AND actual y - 1 is >= 0, read in a scanline to be the fw write
					if (scanlines && (yBlock * yBlockSize + y - wrow) >= 0) {
						band.p_band->RasterIO(
							GF_Read,
							0,
							yBlock * yBlockSize + y - wrow,
							width,
							1,
							p_fwScanline,
							width,
							1,
							band.type,
							0,
							0
						);
					}
					
					size_t blockIndex = scanlines ?
						static_cast<size_t>(y * xBlockSize) :
						static_cast<size_t>(y * xValid);

					//reset no longer used section of focal window matrix
					for (int64_t fwxi = 0; fwxi < fwWidth; fwxi++) {
						focalWindowMatrix[((fwy + wrow - 1) % wrow) * fwWidth + fwxi] = true;
					}

					fwyStart = std::max(fwy, static_cast<int64_t>(0));
					fwyMidStart = std::max(fwy + 1, static_cast<int64_t>(0));
					fwyEnd = scanlines ? 
						std::min(y + yBlock * yBlockSize, height - wrow + 1) :
						std::min(y, yValid - wrow + 1);
					fwyMidEnd = scanlines ? 
						fwyEnd - 1 + static_cast<int64_t>(y + yBlock * yBlockSize > height - wrow + 1) :
						fwyEnd - 1 + static_cast<int64_t>(y > yValid - wrow + 1);
			
					fwx = -wcol + 1;

					for (int x = 0; x < xValid; x++) {
						addSelf = (fwy + verticalPad < 0) ||
							  (y >= yValid - verticalPad) ||
							  (fwx + horizontalPad < 0) ||
							  (x >= xValid - horizontalPad);

						addFw = fwy >= 0 && fwx >= 0;

						if (!scanlines) {
							addSelf &= !(blockLeftPad && x == 0) &&
								   !(blockRightPad && x == xValid - 1) &&
								   !(blockTopPad && y == 0) &&
								   !(blockBotPad && y == yValid - 1);
						}

						fwxStart = std::max(fwx, static_cast<int64_t>(0));
						fwxMidStart = std::max(fwx + 1, static_cast<int64_t>(0));
						fwxEnd = std::min(x, fwWidth);
						fwxMidEnd = scanlines ?
							fwxEnd - 1 + static_cast<int64_t>(y + yBlock * yBlockSize > height - wrow + 1) :
							fwxEnd - 1 + static_cast<int64_t>(y > yValid - wrow + 1);
					
						//GET VAL
						int val = getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex);

						bool isNan = val == nanInt;
						bool accessible = !access.used || (access.used && p_accessMask[blockIndex] == 1);
						
						nextVertSame = !isNan && ((y == yValid - 1) || 
							val == getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex + (scanlines ? xBlockSize : xValid)));
						nextHoriSame = !isNan && ((x == xValid - 1) ||
							val == getPixelValueDependingOnType<int>(band.type, band.p_buffer, blockIndex + 1));

						//UPDATE FOCAL WINDOW MATRIX
						int64_t fwi;
						//set the portion of the focal window matrix to false which is impacted
						//only when the next horizontal pixel value is different than the 
						//current one.
						if (!nextHoriSame) {
							for (int64_t fwyi = fwyStart; fwyi <= fwyMidEnd; fwyi++) {
								fwi = (fwyi % wrow) * fwWidth + fwxEnd;
								focalWindowMatrix[fwi] = false;
							}
						}

						//set the portion of the focal window matrix to false which is impacted
						//only when the next vertical pixel value is different than the current
						//one.
						if (!nextVertSame) {
							for (int64_t fwxi = fwxStart; fwxi <= fwxMidEnd; fwxi++) {
								fwi = (fwyEnd % wrow) * fwWidth + fwxi;
								focalWindowMatrix[fwi] = false;
							}
						}

						//set the portion of the focal window matrix to false which is impacted
						//when either the next vertical or the next horizontal pixel is different
						//than the current one.
						if (!nextHoriSame || !nextVertSame) {
							fwi = (fwyEnd % wrow) * fwWidth + fwxEnd;
							focalWindowMatrix[fwi] = false;
						}

						//set the portion of the focal window matrix to false which must be
						//changed when the next horizontal pixel is different but the previous
						//horizontal pixel is the same as the current one.
						if (!nextHoriSame && prevHoriSame) {
							for (int64_t fwxi = fwxMidStart; fwxi <= fwxMidEnd; fwxi++) {
								fwi = (fwyStart % wrow) * fwWidth + fwxi;
								focalWindowMatrix[fwi] = false;
							}
						}

						//set the portion of the focal window matrix to false which must be
						//changed when the next vertical pixel is different but the previous
						//vertical pixel is the same as the current one. 
						if (!nextVertSame && prevVertSame[x]) {
							for (int64_t fwyi = fwyMidStart; fwyi <= fwyMidEnd; fwyi++) {
								fwi = (fwyi % wrow) * fwWidth + fwxStart;
								focalWindowMatrix[fwi] = false;
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
									focalWindowMatrix[fwi] = false;
								}
							}
						}

						//if the current pixel is nan, or not accessible, mark it's location in 
						//the focal window matrix as false so we know not to add it in the 
						//future without checking for nan or accessibility again.
						if ((isNan || !accessible) && x - horizontalPad >= 0 && y - verticalPad >= 0) {
							int64_t fwyi = ((scanlines ? y + yBlock * yBlockSize : y) - verticalPad) % wrow;
							int64_t fwxi = x - horizontalPad;
							focalWindowMatrix[fwyi * fwWidth + fwxi] = false;
						}

						//ADD INDEX
						if (addSelf && !isNan) {
							strataCounts[val]++;
							if (accessible) {
								Index index = {x + xBlock * xBlockSize, y + yBlock * yBlockSize};

								if (firstXIndexCountPerStrata[val] < firstX) {
									int i = firstXIndexCountPerStrata[val];
									firstXIndexCountPerStrata[val]++;
									firstXIndexesPerStrata[val][i] = index;
								}

								if (randVals[randValIndex]) {
									indexesPerStrata[val].push_back(index);
									indexCountPerStrata[val]++;
								}
								randValIndex++;
							}
						}

						if (addFw) {
							int64_t fwIndex = (fwy % wrow) * fwWidth + fwx;
							
							if (focalWindowMatrix[fwIndex]) {
								queinnecStrataCounts[val]++;
								Index index = {x + xBlock * xBlockSize - wcol, y + yBlock * yBlockSize -wrow};
								int val = scanlines ?
									getPixelValueDependingOnType<int>(band.type, p_fwScanline, x + xBlock * xBlockSize - wcol) :
									getPixelValueDependingOnType<int>(band.type, band.p_buffer, (y - wrow) * xValid + (x - wcol));

								if (queinnecFirstXIndexCountPerStrata[val] < queinnecFirstX) {
									int i = queinnecFirstXIndexCountPerStrata[val];
									queinnecFirstXIndexCountPerStrata[val]++;
									queinnecFirstXIndexesPerStrata[val][i] = index;
								}

								if (queinnecRandVals[queinnecRandValIndex]) {
									queinnecIndexesPerStrata[val].push_back(index);
									queinnecIndexCountPerStrata[val]++;
								}
								queinnecRandValIndex++;
							}
						}

						//INCREMENT WITHIN-BLOCK INDEX
						blockIndex++;
						prevHoriSame = nextHoriSame;
						prevVertSame[x] = nextVertSame;
						fwx++;
					}
					fwy++;
				}	
			}

			//free any firstX vectors which can no longer be used
			for (size_t i = 0; i < numStrata; i++) {
				if (firstXIndexCountPerStrata[i] == firstX) {
					std::vector<Index>().swap(firstXIndexesPerStrata[i]);
					firstXIndexCountPerStrata[i]++;
				}
			}
			
			if (queinnec) {
				for (size_t i = 0; i < numStrata; i++) {
					if (queinnecFirstXIndexCountPerStrata[i] == queinnecFirstX) {
						std::vector<Index>().swap(queinnecFirstXIndexesPerStrata[i]);
						queinnecFirstXIndexCountPerStrata[i]++;
					}
				}
			}

		}
	}

	for (size_t strata = 0; 0 < numStrata; strata++) {
		std::cout << indexesPerStrata.size() << " indexes saved from strata " << strata << std::endl;
		std::cout << "total count: " << strataCounts[strata] << std::endl; 
	}

	if (access.used) {
		if (largeRaster) {
			GDALClose(access.p_dataset);
			free(p_accessMask);
		}
		else {
			free(access.p_dataset);
		}
	}

	uint64_t numDataPixels = 0;
	for (size_t i = 0; i < strataCounts.size(); i++) {
		numDataPixels += strataCounts[i];
	}

	std::vector<uint64_t> strataSampleCounts = calculateAllocation<uint64_t>(
		numSamples,
		allocation,
		strataCounts,
		weights,
		numDataPixels
	);

	std::vector<std::vector<Index> *> strataIndexVectors(numStrata, nullptr);
	std::vector<size_t> nextIndexes(numStrata, 0);

	std::vector<double> xCoords, yCoords;
	std::vector<bool> completedStrata(numStrata, false);
	std::vector<bool> completedStrataQueinnec(numStrata, false);
	std::vector<size_t> samplesAddedPerStrata(numStrata, 0);
	size_t numCompletedStrata = 0;
	size_t numCompletedStrataQueinnec = 0;
	double *GT = p_raster->getGeotransform();
	
	if (queinnec) {
		for (size_t i = 0; i < numStrata; i++) {
			uint64_t totalPixels = queinnecStrataCounts[i];
			uint64_t desiredSamples = strataSampleCounts[i];
			
			uint64_t probIndexesCount = queinnecIndexCountPerStrata[i];
			uint64_t firstXIndexesCount = queinnecFirstXIndexCountPerStrata[i];
			
			if (probIndexesCount >= desiredSamples || firstXIndexesCount < totalPixels) {
				queinnecIndexesPerStrata[i].resize(probIndexesCount);
				auto begin = queinnecIndexesPerStrata[i].begin();
				auto end = queinnecIndexesPerStrata[i].end();
				std::shuffle(begin, end, rng);
				strataIndexVectors[i] =  &queinnecIndexesPerStrata[i];
			}
			else {
				queinnecFirstXIndexesPerStrata[i].resize(firstXIndexesCount);
				auto begin = queinnecFirstXIndexesPerStrata[i].begin();
				auto end = queinnecFirstXIndexesPerStrata[i].end();
				std::shuffle(begin, end, rng);
				strataIndexVectors[i] = &queinnecFirstXIndexesPerStrata[i];
			}
		}

		size_t curStrata = 0;
		while (numCompletedStrataQueinnec < numStrata) {
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

			OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
			p_feature->SetGeometry(&newPoint);
			p_layer->CreateFeature(p_feature);
			OGRFeature::DestroyFeature(p_feature);
	
			if (plot) {
				xCoords.push_back(x);
				yCoords.push_back(y);
			}
	
			curStrata++;
		}
	}
	
	for (size_t i = 0; i < numStrata; i++) {
		uint64_t totalPixels = strataCounts[i];
		uint64_t desiredSamples = strataSampleCounts[i] - samplesAddedPerStrata[i];

		uint64_t probIndexesCount = indexCountPerStrata[i];
		uint64_t firstXIndexesCount = firstXIndexCountPerStrata[i];

		//shuffle the desired vectors so we can iterate starting at 0 without compromising randomness
		if (probIndexesCount >= desiredSamples || firstXIndexesCount < totalPixels) {
			indexesPerStrata[i].resize(probIndexesCount);
			auto begin = indexesPerStrata[i].begin();
			auto end = indexesPerStrata[i].end();
			std::shuffle(begin, end, rng); 
			strataIndexVectors[i] = &indexesPerStrata[i];
		}
		else  {
			firstXIndexesPerStrata[i].resize(firstXIndexesCount);
			auto begin = firstXIndexesPerStrata[i].begin();
			auto end = firstXIndexesPerStrata[i].end();
			std::shuffle(begin, end, rng);
			strataIndexVectors[i] = &firstXIndexesPerStrata[i];
		}

		//set next indexes to 0 if they wre adjusted from adding queinnec indexes
		nextIndexes[i] = 0;
	}
	
	//step 8: generate coordinate points for each sample index.
	size_t curStrata = 0;
	while (numCompletedStrata < numStrata) {
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

		OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
		p_feature->SetGeometry(&newPoint);
		p_layer->CreateFeature(p_feature);
		OGRFeature::DestroyFeature(p_feature);

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
	
/**
 * This function is called by the Python side of the application
 * if the user provided access information. Depending on the
 * method ("random", or "Queinnec"), and depending on the max
 * potential index required, call either strat_queinnec() or 
 * strat_random(), with provided arguments and required template
 * parameters.
 */
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
strat_cpp_access(
	GDALRasterWrapper *p_raster,
	int band,
	size_t numSamples,
	size_t numStrata,
	int wrow,
	int wcol,
	std::string allocation,
	std::string method,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string filename,
	bool largeRaster,
	std::string tempFolder)
{
	return strat(
		p_raster,
		band,
		numSamples,
		numStrata,
		allocation,
		weights,
		method,
		wrow,
		wcol,
		mindist,
		p_access,
		layerName,
		buffInner,
		buffOuter,
		plot,
		filename,
		largeRaster,
		tempFolder
	);
}

/**
 * This function is called by the Python side of the application
 * if the user did not provided access information. Depending on the
 * method ("random", or "Queinnec"), and depending on the max
 * potential index required, call either strat_queinnec() or 
 * strat_random(), with provided arguments and required template
 * parameters.
 */
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
strat_cpp(
	GDALRasterWrapper *p_raster,
	int band,
	size_t numSamples,
	size_t numStrata,
	int wrow,
	int wcol,
	std::string allocation,
	std::string method,
	std::vector<double> weights,
	double mindist,
	bool plot,
	std::string filename,
	bool largeRaster,
	std::string tempFolder)
{
	return strat(
		p_raster, 
		band,
		numSamples,
		numStrata, 
		allocation, 
		weights,
		method,
	       	wrow,			
		wcol,	
		mindist, 
		nullptr, 
		"", 
		0, 
		0,
		plot,	
		filename,
		largeRaster,
		tempFolder
	);
}

PYBIND11_MODULE(strat, m) {
	m.def("strat_cpp", &strat_cpp);
	m.def("strat_cpp_access", &strat_cpp_access);
}
