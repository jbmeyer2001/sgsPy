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
getProbabilityMultiplier(GDALRasterWrapper *p_raster, int numSamples, bool useMindist, double accessibleArea) {
	double height = static_cast<double>(p_raster->getHeight());
	double width = static_cast<double>(p_raster->getWidth());
	double samples = static_cast<double>(numSamples);

	double numer = samples * 2 * (useMindist ? 3 : 1);
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
 * are inaccessable are ignored. The remaining accessable data pixel indexes
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
strat_random(
	GDALRasterWrapper *p_raster,
	int bandNum,
	size_t numSamples,
	size_t numStrata,
	std::string allocation,
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

	std::vector<unsigned long long int> strataCounts(numStrata, 0);
	std::vector<std::vector<Index>> indexesPerStrata(numStrata);
	std::vector<size_t> indexCountPerStrata(numStrata, 0);

	size_t firstX = 10000;
	std::vector<std::vector<Index>> firstXIndexesPerStrata(numStrata, std::vector<Index>(firstX));
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
	//and accessable area, we can estimate a percentage chance for each pixel and set up a multiplier
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
	uint64_t multiplier = getProbabilityMultiplier(p_raster, numSamples, useMindist, access.area);
	
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;

			//READ BLOCK
			//bandMutex.lock();
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			rasterBandIO(band, band.p_buffer, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true); 
			//bandMutex.unlock();									//read = true

			//READ ACCESS BLOCK IF NECESSARY
			if (access.used) {
				//accessMutex.lock();
				rasterBandIO(access.band, p_accessMask, xBlockSize, yBlockSize, xBlock, yBlock, xValid, yValid, true);
				//accessMutex.lock();									      //read = true	
			}
			
			//CALCULATE RAND VALUES
			//rngMutex.lock();
			for (int i = 0; i < randValIndex; i++) {
				randVals[i] = ((rng() >> 11) & multiplier) == multiplier;
			}
			randValIndex = 0;
			//rngMutex.unlock();

			//ITERATE THROUGH BLOCK AND UPDATE VECTORS
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

					//std::cout << "HERE 4.7" << std::endl;
					//std::cout << "index.x: " << index.x << std::endl;
					//std::cout << "index.y: " << index.y << std::endl;
					//std::cout << "val: " << val << std::endl;
					//std::cout << "strataCounts.size(): " << strataCounts.size() << std::endl;
					//std::cout << "indexesPerStrata.size(): " << indexesPerStrata.size() << std::endl;
					//std::cout << "indexesPerStrata[val].size(): " << indexesPerStrata[val].size() << std::endl;
					//std::cout << "indexCountPerStrata[val]: " << indexCountPerStrata[val] << std::endl;

					//UPDATE VECTORS
					strataCounts[val]++;
					indexesPerStrata[val].push_back(index);
					indexCountPerStrata[val]++;

					//std::cout << "HERE 4.8" << std::endl;

					//std::cout << "HERE 4.9" << std::endl;
					//INCREMENT WITHIN-BLOCK INDEX
					blockIndex++;
				}
			}


			//free any firstX vectors which can no longer be used
			for (size_t i = 0; i < numStrata; i++) {
				if (firstXIndexCountPerStrata[i] == firstX) {
					std::vector<Index>().swap(firstXIndexesPerStrata[i]);
					firstXIndexCountPerStrata[i]++;
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
	for (const uint64_t& count : strataCounts) {
		numDataPixels += count;
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
	for (size_t i = 0; i < numStrata; i++) {
		uint64_t totalPixels = strataCounts[i];
		uint64_t desiredSamples = strataSampleCounts[i];

		uint64_t probIndexesCount = indexCountPerStrata[i];
		uint64_t firstXIndexesCount = firstXIndexCountPerStrata[i];

		//shuffle the desired vectors so we can iterate starting at 0 without compromising randomness
		if (probIndexesCount >= desiredSamples || firstXIndexesCount < totalPixels) {
			auto begin = indexesPerStrata[i].begin();
			auto end = indexesPerStrata[i].end();
			std::shuffle(begin, end, rng); 
			strataIndexVectors[i] = &indexesPerStrata[i];
		}
		else  {
			auto begin = firstXIndexesPerStrata[i].begin();
			auto end = firstXIndexesPerStrata[i].end();
			std::shuffle(begin, end, rng);
			strataIndexVectors[i] = &firstXIndexesPerStrata[i];
		}
	}
	
	//step 8: generate coordinate points for each sample index.
	std::vector<double> xCoords, yCoords;
	std::vector<bool> completedStrata(numStrata, false);
	std::vector<size_t> samplesAddedPerStrata(numStrata, 0);
	size_t numCompletedStrata = 0;
	double *GT = p_raster->getGeotransform();
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
 * This function conducts stratified sampling on the provided stratified raster
 * using the Queinnec method. The queinnec method first tries to sample pixels which 
 * are surrounded entirely by pixels of the same strata within a user-defined 
 * focal window. If there are not enough focal-window pixels to sample the desired number,
 * pixels are randomly sampled.
 *
 *
 * First, the raster is iterated over, the no data pixels and the pixels which
 * are inaccessable are not considered for sampling, however the inaccessable
 * pixels are still considered as part of a valid focal window for accessable pixels.
 *
 * The remaining accessable data pixel indexes are placed into either 
 * queinnec-sampling vectors or random-sampling vectors corresponding to their strata. 
 *
 * Rather than checking the surrounding focal window pixels for every single pixel --
 * a computationally expensive method -- a moving focal window of boolean values
 * (the 'focal window matrix') is used to determine whether a pixel should go into
 * the queinnec vectors or the random vectors. The focal window matrix is initialized
 * to 'true', and as the stratified raster is iterated though, each pixel checks the
 * next vertical and next horizontal pixel to see if they are the same -- boolean 
 * values for the previous pixels are also saved. Then, using these 4 booleans (
 * 'next horizontal same', 'next vertical same', 'previous horizontal same',
 * and 'previous vertical same') the corresponding pixels in the focal window
 * matrix are set to false if required. When all of the pixels which may
 * influence a given focal window pixel are checked, that pixel is added to
 * the queinnec vectors if the focal window value is true, and the
 * random vectors otherwise.
 *
 * Due to the relatively complex nature of the focal window matrix, inline
 * tests are included to ensure the correct pixels are added to the queinnec
 * vectors and random vectors. They are commented out, but are left there
 * for if the implementation requires future adjustments.
 *
 *
 * Next, the calculate_allocation function is used to determine the the total
 * number of samples which should be allocated to each strata, depending on the
 * number of pixels in each strata and the allocation method specified by the
 * user.
 *
 *
 * Finally, random indexes in the queinnec strata vectors are selected, 
 * and using their index value in the stratified raster and geotransform, 
 * a point geometry is calculated.
 *
 * When a stratum is allocated more than half the number of pixels in the 
 * queinnec strata vector, random indexes in the queinnec strata vectors
 * are selected which WONT be included, and the remaining are iterated over
 * in a random order.
 *
 * When mindist is inequal to zero, three times the number of samples are 
 * randomly selected in order to ensure a number of samples as close to the 
 * desired number are selected, as having a large mindist may mean some samples
 * can't be included due to their proximity to other already-selected pixels.
 *
 * If a strata does not have the desired number of samples, the random
 * strata vectors are then used to added any required remaining
 * samples using the same random selection technique. 
 */

std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
strat_queinnec(
	GDALRasterWrapper *p_raster,
	int band,
	size_t numSamples,
	size_t numStrata,
	int wrow,
	int wcol,
	std::string allocation,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string filename)
{
	//step 1: get raster band
	GDALDataType type = p_raster->getRasterBandType(band);
	printTypeWarningsForInt32Conversion(type);	

	GDALRasterBand *p_band = p_raster->getRasterBand(band);
	void *p_data = VSIMalloc3(
		p_raster->getHeight(),
		p_raster->getWidth(),
		p_raster->getRasterBandTypeSize(band)
	);
	CPLErr err = p_band->RasterIO(
		GF_Read,
		0,
		0,
		p_raster->getWidth(),
		p_raster->getHeight(),
		p_data,
		p_raster->getWidth(),
		p_raster->getHeight(),
		type,
		0,
		0
	);
	if (err) {
		throw std::runtime_error("error reading raster band from dataset.");
	}
	
	//step 2: create stratum index storing vectors
	std::vector<std::vector<size_t>> queinnecStratumIndexes;
	std::vector<std::vector<size_t>> randomStratumIndexes;
	queinnecStratumIndexes.resize(numStrata);
	randomStratumIndexes.resize(numStrata);

	//step 3: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	uint8_t *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first;
		p_mask = static_cast<uint8_t *>(maskInfo.second);
	}
	
	//step 4: determine size of and allocate focal window matrix
	size_t width = p_raster->getWidth();
	size_t height = p_raster->getHeight();
	size_t horizontalPad = wcol / 2;
	size_t verticalPad = wrow / 2;
	size_t fwHeight = wrow;
	size_t fwWidth = width - wcol + 1;
	std::vector<bool> focalWindowMatrix(fwWidth * fwHeight, true);
	std::vector<bool> prevVertSame(width, true);
	bool prevHoriSame = true;
	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
	size_t noDataPixelCount = 0;
	bool nextVertSame;
	bool nextHoriSame;

	size_t y = 0;
	int64_t fwy = -wrow + 1;
	size_t x;
	int64_t fwx;
	bool addSelf;
	bool addfw;

	//step 5: iterate through strat raster.
	while (y < height) {
		//reset no longer used section of focal window matrix
		for (int64_t fwxi = 0; fwxi < fwWidth; fwxi++) {
			focalWindowMatrix[((fwy + wrow - 1) % wrow) * fwWidth + fwxi] = true;
		}

		int64_t fwyStart = std::max(fwy, static_cast<int64_t>(0));
		int64_t fwyMidStart = std::max(fwy + 1, static_cast<int64_t>(0));
		int64_t fwyEnd = std::min(y, static_cast<size_t>(height - wrow + 1));

		int64_t fwyMidEnd = fwyEnd - 1 + static_cast<int64_t>(y > static_cast<size_t>(height - wrow + 1));

		x = 0;
		fwx = -wcol + 1;

		while (x < width) {
			addSelf = (fwy + verticalPad < 0) || 
				  (y >= height - verticalPad) || 
				  (fwx + horizontalPad < 0) || 
				  (x >= width - horizontalPad);
			addfw = fwy >= 0 && fwx >= 0;

			int64_t fwxStart = std::max(fwx, static_cast<int64_t>(0));
			int64_t fwxMidStart = std::max(fwx + 1, static_cast<int64_t>(0));
			int64_t fwxEnd = std::min(x, fwWidth);
			int64_t fwxMidEnd = fwxEnd - 1 + static_cast<int64_t>(y > static_cast<size_t>(height - wrow + 1));

			size_t index = y * width + x;
			int val = getPixelValueDependingOnType<int>(type, p_data, index);
			bool isNan = std::isnan(val) || static_cast<double>(val) == noDataValue;
			bool accessable = !p_access || (p_mask[index] != 0); //check access mask
			noDataPixelCount += (size_t)(isNan || !accessable);
			
			//allow inaccessible pixels to count towards focal window
			nextVertSame = !isNan && ((y == height - 1) || 
					val == getPixelValueDependingOnType<int>(type, p_data, index + width));
		       	nextHoriSame = !isNan && ((x == width - 1) || 
					val == getPixelValueDependingOnType<int>(type, p_data, index + 1));	

			//add the current pixel if it is not within the focal window matrix (too close to raster edges)
			if (addSelf && !isNan && accessable) {
				randomStratumIndexes[(size_t)val].push_back(index);
			}

			//add the focal window pixel which can no longer be altered by future pixel values
			if (addfw) {
				int64_t fwIndex = (fwy % wrow) * fwWidth + fwx;
				index = (fwy + verticalPad) * width + fwx + horizontalPad;	
				val = getPixelValueDependingOnType<int>(type, p_data, index);
				isNan = std::isnan(val) || static_cast<double>(val) == noDataValue;
				accessable = (!p_access) || (p_mask[index] != 0); //check access mask

				if (!isNan && accessable && focalWindowMatrix[fwIndex]) {
					queinnecStratumIndexes[(size_t)val].push_back(index);
				}
				else if (!isNan && accessable) {
					randomStratumIndexes[(size_t)val].push_back(index);
				}

			}

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
			
			prevHoriSame = nextHoriSame;
			prevVertSame[x] = nextVertSame;
			x++;
			fwx++;
		}

		y++;
		fwy++;
	}
	if (p_access) {
		free(p_accessMaskDataset);	
	}

	size_t numDataPixels = p_raster->getWidth() * p_raster->getHeight() - noDataPixelCount;

	//step 6: test to ensure focal queinnec/random vector additions are correct
	/*
	//FOR TESTING PURPOSES	
	size_t checkNumDataPixels = 0;
	for (size_t i = 0; i < numStrata; i++) {
		checkNumDataPixels += queinnecStratumIndexes[i].size() + randomStratumIndexes[i].size();
	}

	if (checkNumDataPixels != numDataPixels) {
		throw std::runtime_error("**BUG** incorrect number of pixels added to random and queinnec stratum indexes.");
	}

	//FOR TESTING PURPOSES
	for (size_t i = 0; i < queinnecStratumIndexes.size(); i++) {
		for (size_t j = 0; j < queinnecStratumIndexes[i].size(); j++) {
			U index = queinnecStratumIndexes[i][j];
			for (int ii = 0; ii < wrow; ii++) {
				for (int jj = 0; jj < wcol; jj++) {
					U checkIndex = index + (ii - (wrow / 2)) * width + (jj - (wcol / 2));
					if (p_strata[index] != p_strata[checkIndex]) {
						throw std::runtime_error("**BUG** index " + std::to_string(checkIndex) 
							+ " different than " + std::to_string(index) + ".");
					}
				}
			}
		}
	}

	//FOR TESTING PURPOSES
	for (size_t i = 0; i < randomStratumIndexes.size(); i++) {
		for (size_t j = 0; j < randomStratumIndexes[i].size(); j++) {
			int64_t index =  randomStratumIndexes[i][j];
			bool same = true;
			for (int ii = 0; ii < wrow; ii++) {
				for (int jj = 0; jj < wcol; jj++) {
					int64_t checkIndex = index + (ii - (wrow / 2)) * width + (jj - (wcol / 2));
					if (checkIndex >= 0 && checkIndex < width * height) {
						same &= (p_strata[index] != p_strata[checkIndex]);
					}
					else {
						same = false;
					}
				}
			}

			if (same) {
				throw std::runtime_error("**BUG** index " + std::to_string(index) + " is surrounded by all same pixels.");
			}
		}
	}
	*/

	//step 7: calculate allocation of samples depending on stratum sizes 
	std::vector<size_t> strataSizes;
	for (size_t i = 0; i < queinnecStratumIndexes.size(); i++) {
		strataSizes.push_back(randomStratumIndexes[i].size() + queinnecStratumIndexes[i].size());
	}

	std::vector<size_t> stratumCounts = calculateAllocation<size_t>(
		numSamples,
		allocation,
		strataSizes,
		weights,
		numDataPixels
	);

	//step 8: create new in-memory dataset to store sample points
	GDALAllRegister();
	GDALDataset *p_sampleDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	OGRLayer *p_sampleLayer = p_sampleDataset->CreateLayer("samples", nullptr, wkbPoint, nullptr);

	//step 9: determine queinnec indexes to try including as samples
	std::vector<std::unordered_set<size_t>> sampleIndexes;
	std::vector<typename std::unordered_set<size_t>::iterator> sampleIterators;
	std::vector<size_t> samplesAdded;
	std::vector<size_t> strataNum;
	size_t totalSamplesAdded = 0;

	std::mt19937::result_type seed = time(nullptr);
	for (size_t i = 0; i < stratumCounts.size(); i++) {	
		sampleIndexes.push_back({});
		samplesAdded.push_back(0);
		strataNum.push_back(i);

		auto rng = std::bind(
			std::uniform_int_distribution<size_t>(0, queinnecStratumIndexes[i].size() - 1),
			std::mt19937(seed)
		);

		size_t stratumSamples = std::min((mindist == 0) ? (size_t)stratumCounts[i] : (size_t)stratumCounts[i] * 3, queinnecStratumIndexes[i].size());

		if (stratumSamples > queinnecStratumIndexes[i].size() / 2) {
			std::unordered_set<size_t> dontSamplePixels;
			while (dontSamplePixels.size() < queinnecStratumIndexes[i].size() - stratumSamples) {
				dontSamplePixels.insert(rng());
			}
			for (size_t j = 0; j < queinnecStratumIndexes[i].size(); j++) {
				if (dontSamplePixels.find(j) == dontSamplePixels.end()) {
					sampleIndexes[i].insert(queinnecStratumIndexes[i][j]);
				}
			}
		}
		else {
			while (sampleIndexes[i].size() < stratumSamples) {
				sampleIndexes[i].insert(queinnecStratumIndexes[i][rng()]);
			}
		}

		sampleIterators.push_back(sampleIndexes[i].begin());
	}

	//step 10: generate coordinate points for each queinnec sample, and only add if they're outside of mindist
	std::vector<double> xCoords, yCoords;

	size_t sIndex = 0;
	size_t completedStratum = 0;
	double *GT = p_raster->getGeotransform();
	while (completedStratum < stratumCounts.size()) {
		if (sIndex >= strataNum.size()) {
			sIndex = 0;
		}
		size_t strata = strataNum[sIndex];

		if (sampleIterators[strata] == sampleIndexes[strata].end() || stratumCounts[strata] == samplesAdded[strata]) {
			completedStratum++;
			strataNum.erase(strataNum.begin() + sIndex);
		}
		else {
			size_t index = *sampleIterators[strata];
			sampleIterators[strata] = std::next(sampleIterators[strata]);

			double yIndex = index / p_raster->getWidth();
			double xIndex = index - (yIndex * p_raster->getWidth());
			double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
			double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
			OGRPoint newPoint = OGRPoint(xCoord, yCoord);
	
			if (mindist != 0.0 && p_sampleLayer->GetFeatureCount() != 0) {
				bool add = true;
				for (const auto& p_feature : *p_sampleLayer) {
					OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
					if (newPoint.Distance(p_point) < mindist) {
						add = false;
						break;
					}
				}
				if (!add) {
					sIndex++;
					continue;
				}
			}

			OGRFeature *p_feature = OGRFeature::CreateFeature(p_sampleLayer->GetLayerDefn());
			p_feature->SetGeometry(&newPoint);
			p_sampleLayer->CreateFeature(p_feature);
			OGRFeature::DestroyFeature(p_feature);

			samplesAdded[strata]++;
			totalSamplesAdded++;
			sIndex++;

			if (plot) {
				xCoords.push_back(xCoord);
				yCoords.push_back(yCoord);
			}
		}
	}

	//step 11: determine random indexes to try including as samples
	for (size_t i = 0; i < stratumCounts.size(); i++) {	
		stratumCounts[i] -= samplesAdded[i];
		if (stratumCounts[i] > 0) {
			completedStratum--;
			sampleIndexes[i] = {};
			samplesAdded[i] = 0;
			strataNum.push_back(i);

			auto rng = std::bind(
				std::uniform_int_distribution<size_t>(0, randomStratumIndexes[i].size() - 1),
				std::mt19937(seed)
			);

			size_t stratumSamples = std::min((mindist == 0) ? (size_t)stratumCounts[i] : (size_t)stratumCounts[i] * 3, randomStratumIndexes[i].size());

			if (stratumSamples > randomStratumIndexes[i].size() / 2) {
				std::unordered_set<size_t> dontSamplePixels;
				while (dontSamplePixels.size() < randomStratumIndexes[i].size() - stratumSamples) {
					dontSamplePixels.insert(rng());
				}
				for (size_t j = 0; j < randomStratumIndexes[i].size(); j++) {
					if (dontSamplePixels.find(j) == dontSamplePixels.end()) {
						sampleIndexes[i].insert(randomStratumIndexes[i][j]);
					}
				}
			}
			else {
				while (sampleIndexes[i].size() < stratumSamples) {
					sampleIndexes[i].insert(randomStratumIndexes[i][rng()]);
				}
			}

			sampleIterators[i] = sampleIndexes[i].begin();
		}
	}
	//step 12: try adding random samples
	sIndex = 0;
	while (completedStratum < stratumCounts.size()) {
		//determine strata from i
		if (sIndex >= strataNum.size()) {
			sIndex = 0;
		}
		size_t strata = strataNum[sIndex];

		if (sampleIterators[strata] == sampleIndexes[strata].end() || stratumCounts[strata] == samplesAdded[strata]) {
			completedStratum++;
			strataNum.erase(strataNum.begin() + sIndex);
		}
		else {
			size_t index = *sampleIterators[strata];
			sampleIterators[strata] = std::next(sampleIterators[strata]);

			double yIndex = index / p_raster->getWidth();
			double xIndex = index - (yIndex * p_raster->getWidth());
			double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
			double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
			OGRPoint newPoint = OGRPoint(xCoord, yCoord);
	
			if (mindist != 0.0 && p_sampleLayer->GetFeatureCount() != 0) {
				bool add = true;
				for (const auto& p_feature : *p_sampleLayer) {
					OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
					if (newPoint.Distance(p_point) < mindist) {
						add = false;
						break;
					}
				}
				if (!add) {
					sIndex++;
					continue;
				}
			}

			OGRFeature *p_feature = OGRFeature::CreateFeature(p_sampleLayer->GetLayerDefn());
			p_feature->SetGeometry(&newPoint);
			p_sampleLayer->CreateFeature(p_feature);
			OGRFeature::DestroyFeature(p_feature);

			samplesAdded[strata]++;
			totalSamplesAdded++;
			sIndex++;

			if (plot) {
				xCoords.push_back(xCoord);
				yCoords.push_back(yCoord);
			}
		}
	}

	//step 13: create GDALVectorWrapper with dataset containing points
	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_sampleDataset);

	//step 14: write to file if filename is not ""
	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	return {{xCoords, yCoords}, p_sampleVectorWrapper, totalSamplesAdded};
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
	if (method == "random") {
		return strat_random(
			p_raster,
			band,
			numSamples,
			numStrata,
			allocation,
			weights,
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
	else { //method == "Queinnec"
		return strat_queinnec(
			p_raster,
			band,
			numSamples,
			numStrata,
			wrow,
			wcol,
			allocation,
			weights,
			mindist,
			p_access,
			layerName,
			buffInner,
			buffOuter,
			plot,
			filename
		);
	}
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
	if (method == "random") {
		return strat_random(
			p_raster, 
			band,
			numSamples,
			numStrata, 
			allocation, 
			weights, 
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
	else { //method == "Queinnec"
		return strat_queinnec(
			p_raster,
			band,
			numSamples,
			numStrata,
			wrow,
			wcol,
			allocation,
			weights,
			mindist,
			nullptr,
			"",
			0,
			0,
			plot,
			filename
		);
	}
}

PYBIND11_MODULE(strat, m) {
	m.def("strat_cpp", &strat_cpp);
	m.def("strat_cpp_access", &strat_cpp_access);
}
