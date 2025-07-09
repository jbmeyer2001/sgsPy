/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of stratified sampling
 * Author: Joseph Meyer
 * Date: June, 2025

 *
 ******************************************************************************/

#include <iostream>
#include <random>

#include "access.h"
#include "raster.h"
#include "vector.h"
#include "write.h"

/**
 * Helper function which calculates the count of values for each strata.
 *
 * @param std::string allocation method
 * @param std::vector<size_t>& sizes of each strata
 * @param std::vector<double> weights of each strata
 * @param size_t total number of pixels
 * @returns std::vector<size_t> counts of each stratum
 */
template <typename U>
std::vector<U>
calculateAllocation(
	size_t numSamples,
	std::string allocation, 
	std::vector<U> *p_strataSizes, 
	std::vector<double> *p_weights,
	U numPixels)
{
	std::vector<U> retval;
	size_t remainder;
	size_t numStrata = p_strataSizes->size();
	if (allocation == "prop") {
		//allocate the samples per stratum according to stratum size
		remainder = numSamples;
		for (size_t i = 0; i < numStrata; i++) {
			size_t count = numSamples / (numPixels / p_strataSizes->at(i));
			retval.push_back(count);
			remainder -= count;
		}

		//redistribute remainder pixels equally among strata
		size_t i = 0;
		while (remainder > 0) {
			retval[i] += 1;
			remainder -= 1;
			i++;
		}
	}
	else if (allocation == "equal") { //TODO what if a strata doesn't have enough pixels???
		//determine the count of samples per strata
		size_t strataSampleCount = numSamples / numStrata;
		
		for (size_t i = 0; i < numStrata; i++) {
			if (p_strataSizes->at(i) < strataSampleCount) {
				retval.push_back(p_strataSizes->at(i));
				std::cout << "warning: strata " << i << " does not have enough pixels for the full " << strataSampleCount << " samples it should recieve. There will be less than " << numSamples << " final samples." << std::endl;
			}
			else {
				retval.push_back(strataSampleCount);
			}

		}
	}
	else if (allocation == "manual") { //TODO what if a strata doesn't have enough pixels???
		if (numStrata != p_weights->size()) {
			throw std::runtime_error("size of weights is not equal to the number of strata.");
		}

		//allocate samples accordign to weights.
		remainder = numSamples;
		for (size_t i = 0; i < numStrata; i++) {
			size_t count = (size_t)((double)p_strataSizes->at(i) * p_weights->at(i));
			retval.push_back(count);
			remainder -= count;
		}

		//redistribute remainder pixels equally among strata
		size_t i = 0;
		while (remainder > 0) {
			retval[i] += 1;
			remainder -= 1;
			i++;
		}

		//ensure number of strata doesn't exceed maximum
		for (size_t i = 0; i < retval.size(); i++) {
			if (retval[i] > p_strataSizes->at(i)) {
				std::cout << "warning: strata " << i << " does not have enough pixels for the full " << retval[i] << " samples it should recieve. There will be less than " << numSamples << " final samples." << std::endl;
				retval[i] = p_strataSizes->at(i);
			}
		}
	}
	else { //allocaiton == "optim"
		//TODO implement
		throw std::runtime_error("'optim' has not been implemented!");
	}	

	return retval;
}

template <typename U>
inline bool checkNan(U index, U& val, uint8_t *p_mask, float *p_strat, double noDataValue) {
	bool isNan = p_mask && p_mask[index] == 0;
	if (!isNan) {
		val = p_strat[index];
		isNan = std::isnan(val) || (double)val = noDataValue;
	}
	return isNan;
}

/**
 *
 */
template <typename U>
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_random(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	U numStrata,
	std::string allocation,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string filename)
{
	//step 1: get raster band
	if (p_raster->getRasterType() != GDT_Float32) {
		throw std::runtime_error("raster band must be of type float32.");
	}
	float *p_strat = (float *)p_raster->getRasterBand(0);

	//step 2: create stratum index storing vectors
	std::vector<std::vector<U>> stratumIndexes;
	stratumIndexes.resize(numStrata);

	//step 3: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using mask
		p_mask = maskInfo.second; //pointer to mask
	}

	//stpe 4: iterate through raster band
	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
	U noDataPixelCount = 0;
	for (U i = 0; i < (U)p_raster->getWidth() * (U)p_raster->getHeight(); i++) {
		if (p_access) {
			//step 4.1 if the current pixel is not accessable, mark it as nodata and don't read it
			if (((uint8_t *)p_mask)[i] == 0) {
				noDataPixelCount++;
				continue;
			}
		}

		float val = p_strat[i];
		if (std::isnan(val) || (double)val == noDataValue) {
			//step 4.2: increment noDataPixelCount if encountered noData
			noDataPixelCount++;
		}
		else {
			//step 4.3: add data index to stratum index array
			stratumIndexes[(size_t)val].push_back(i);
		}
	}
	if (p_access) {
		free(p_accessMaskDataset);
	}

	U numDataPixels = (p_raster->getWidth() * p_raster->getHeight()) - noDataPixelCount;

	//step 5: using the size of the stratum and allocation method,
	//determine the number of samples to take from each stratum
	std::vector<U> strataSizes;
	for (size_t i = 0; i < stratumIndexes.size(); i++) {
		strataSizes.push_back(stratumIndexes[i].size());
	}

	std::vector<U> stratumCounts = calculateAllocation<U>(
		numSamples,
		allocation,
		&strataSizes,
		&weights,
		numDataPixels
	);
	
	//Step 6: generate random number generator using mt19937	
	std::vector<std::function<U(void)>> rngs;
	for (size_t i = 0; i < stratumCounts.size(); i++) {
		std::mt19937::result_type seed = time(nullptr);
		rngs.push_back(std::bind(
			std::uniform_int_distribution<U>(0, stratumCounts[i] - 1),
			std::mt19937(seed)
		));
	}

	//step 7: determine pixel values to include as samples.
	
	//make a set of pixels to include as samples
	std::unordered_set<U> samplePixels;

	//make backup unordered set of pixels to use (if mindist causes others not to be included);
	std::unordered_set<U> backupSamplePixels;

	for (size_t i = 0; i < stratumCounts.size(); i++) {
		size_t startSize = samplePixels.size();
		while(samplePixels.size() < startSize + stratumCounts[i]) {
			samplePixels.insert(stratumIndexes[i][rngs[i]()]);
		}
		if (mindist != 0.0) {
			size_t startSize = backupSamplePixels.size();
			//TODO what if there aren't enough pixels???
			while(backupSamplePixels.size() < startSize + (stratumCounts[i] * 2)) {
				U pixel = stratumIndexes[i][rngs[i]()];
				if (samplePixels.find(pixel) == samplePixels.end()) {
					backupSamplePixels.insert(pixel);
				}
			}
		}
	}

	//step 8: generate coordinate points for each sample index, and only add if they're outside of mindist
	std::vector<double> xCoords;
	std::vector<double> yCoords;
	std::vector<OGRPoint> points;
	std::vector<std::string> wktPoints;

	double *GT = p_raster->getGeotransform();
	for (auto index : samplePixels) {
		//TODO check if we can move this to a helper function...
		double yIndex = index / p_raster->getWidth();
		double xIndex = index - (yIndex * p_raster->getWidth());
		double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
		double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
		OGRPoint newPoint = OGRPoint(xCoord, yCoord);
	
		if (mindist != 0.0 && points.size() != 0) {
			U pIndex = 0;
			while ((pIndex < points.size()) && (newPoint.Distance(&points[pIndex]) > mindist)) {
				pIndex++;
			}
			if (pIndex != (U)points.size()) {
				continue;
			}
		}
		points.push_back(newPoint);
		wktPoints.push_back(newPoint.exportToWkt());
		xCoords.push_back(xCoord);
		yCoords.push_back(yCoord);
	}
	
	if (mindist != 0.0 && points.size() < samplePixels.size()) {
		for (auto index : backupSamplePixels) {
			double yIndex = index / p_raster->getWidth();
			double xIndex = index - (yIndex * p_raster->getWidth());
			double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
			double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
			OGRPoint newPoint = OGRPoint(xCoord, yCoord);

			U pIndex = 0;
			while ((pIndex < points.size()) && (newPoint.Distance(&points[pIndex]) > mindist)) {
				pIndex++;
			}

			if (pIndex == points.size()) {
				points.push_back(newPoint);
				wktPoints.push_back(newPoint.exportToWkt());
				xCoords.push_back(xCoord);
				yCoords.push_back(yCoord);
			}

			if (points.size() == numSamples) { 
				break;
			}
		}
	}

	if (filename != "") {
		try {
			writeSamplePoints(points, filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	return {{xCoords, yCoords}, wktPoints};
}
	
/**
 *
 */
template <typename U>
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_queinnec(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	U numStrata,
	int wrow,
	int wcol,
	std::string allocation,
	std::vector<double> weights,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string fileName)
{
	//step 1: get raster band
	if (p_raster->getRasterType() != GDTFloat32) {
		throw std::runtime_error("raster band must be of type float32.");
	}
	float *p_strat = p_raster->getBand(0);

	//step 2: create stratum index storing vectors
	std::vector<std::vector<U>> queinnecStratumIndexes;
	std::vector<std::vector<U>> stratumIndexes;
	queinnecStratumIndexes.resize(numStrata);
	randomStratumIndexes.resize(numStrata);

	//step 3: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first;
		p_mask = maskInfo.second;
	}	

	int width = p_raster->getWidth();
	int height = p_raster->getHeight();
	int horizontalPad = (wcol / 2);
	int verticalPad = (wrow / 2);
	U fwMatrixHeight = wrow;
	U fwMatrixWidth = width - wcol + 1;

	std::vector<bool> focalWindowMatrix(true, fwMatrixWidth * fwMatrixHeight);
	std::vector<bool> prevVertSame(true, width);
	bool prevHoriSame = true;

	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
	U noDataPixelCount = 0;

	bool addSelf;
	bool addUpperLeftCornerPixel;

	U fwUpperLeftX;
	U fwUpperLeftY = 0 - verticalPad * 2;

	U fwMatrixXStart;
	U fwMatrixXEnd;
	U fwMatrixYStart;
	U fwMatrixYEnd = 0;
		
	for (U y = 0; y < height; y++) {
		addSelf = y < verticalPad || y > height - 1 - verticalPad;
		addUpperLeftCornerPixel = y >= verticalPad * 2;
		
		//fwUpperLeftY and fwMatrixYEnd already calculated/set elsewhere
		fwMatrixYStart = std::max(fwUpperLeftY + 1, 0);
		if (fwMatrixYStart == fwMatrixHeight) {
			fwMatrixYStart = 0;
		}

		for (x = 0; x < width; x++) {
			addSelf |= x < horizontalPad || x > width - 1 - horizontalPad;
			addUpperLeftCornerPixel &= x >= horizontalPad * 2;

			fwUpperLeftX = x - horizontalPad * 2;
			fwMatrixXEnd = std::min(x, fwMatrixWidth - 1);
			fwMatrixXStart = std::max(0, fwUpperLeft + 1);

			U index = y * width + x;
			U val;
			bool isNan = checkNan(index, &val, p_mask, p_strata, noDataValue);
		
			if (addSelf && !isNan) {
				randomStratumIndexes[(size_t)val].push_back(index);			
			}

			if (addUpperLeftCornerPixel) {
				U upperLeftIndex = (y - verticalPad) * width + (x - horizontalPad);
				U upperLeftfwIndex = fwUpperLeftY * fwMatrixWidth + fwUpperLeftX;
				U upperLeftVal;
				bool upperLeftIsNan = checkNan(upperLeftIndex, &upperLeftVal, p_mask, p_strata, noDataValue);

				//if the corrosponding value in the focal window matrix is true, it means
				//that every surrounding pixel is the same, so the pixel should be added
				//to the queinnec group of pixels, otherwise add to the random group of pixels
				if (focalWindowMatrix[upperLeftfwIndex] && !upperLeftIsNan) {
					queinnecStratumIndexes[(size_t)upperLeftVal].push_back(upperLeftIndex);
				}
				else if (!upperLeftIsNan) {
					randomStratumIndexes[(size_t)upperLeftVal].push_back(upperLeftIndex);
				}
			}

			nextVertSame == !isNan && ((y == height - 1) || val == p_strata[index + width]);
			nextHoriSame == !isNan && ((x == width - 1) || val == p_strata[index + 1]);
			
			//add special case for if previous vertical and horizontal were the same
			//but next vertical and horizontal were different, to avoid rewriting
			//focal window matrix pixels.
			if (prevVertSame[x] && prevHoriSame && !nextVertSame && !nextHoriSame) {
				//TODO add
				prevVertSame = false;
				prevHoriSame = false;
				continue;
			}

			//conceptually split the affected pixels in the focal window matrix into
			//pixels which need to be updated if:
			//
			//either !nextHoriSame OR !nextVertSame
			//!nextHoriSame
			//!nextVertSame
			//

			prevVertSame[x] = nextVertSame;
			prevHoriSame = nextHoriSame;
		}

		//TODO set required values in focal window matrix to true
		fwUpperLeftY++;
		if (fwUpperLeftY == fwMatrixHeight) {
			fwUpperLeftY -= 0;
		}

		fwMatrixYEnd++;
		if (fwMatrixYEnd == fwMatrixHeight) {
			fwMatrixYEnd++;
		}
	}
			
	return {{{0.0}, {0.0}}, {""}};
}

/**
 *
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_cpp_access(
	GDALRasterWrapper *p_raster,
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
	std::string filename)
{
	std::string minIndexIntType = p_raster->getMinIndexIntType(true); //single band
	if (method == "random") {
		if (minIndexIntType == "unsigned_short") { return strat_random<unsigned short>(p_raster,numSamples,numStrata,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
		if (minIndexIntType == "unsigned") { return strat_random<unsigned>(p_raster,numSamples,numStrata,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
		if (minIndexIntType == "unsigned_long") { return strat_random<unsigned long>(p_raster,numSamples,numStrata,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
		if (minIndexIntType == "unsigned_long_long") { return strat_random<unsigned long long>(p_raster,numSamples,numStrata,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
	}
	else { //method == "Queinnec"
		if (minIndexIntType == "unsigned_short") { return strat_queinnec<unsigned short>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
		if (minIndexIntType == "unsigned") { return strat_queinnec<unsigned>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
		if (minIndexIntType == "unsigned_long") { return strat_queinnec<unsigned long>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
		if (minIndexIntType == "unsigned_long_long") { return strat_queinnec<unsigned long long>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,p_access,layerName,buffInner,buffOuter,filename); }
	}
}

/**
 *
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_cpp(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	size_t numStrata,
	int wrow,
	int wcol,
	std::string allocation,
	std::string method,
	std::vector<double> weights,
	double mindist,
	std::string filename)
{
	std::string minIndexIntType = p_raster->getMinIndexIntType(true); //single band
	if (method == "random") {
		if (minIndexIntType == "unsigned_short") { return strat_random<unsigned short>(p_raster,numSamples,numStrata,allocation,weights,mindist,nullptr,"",0,0,filename); }
		if (minIndexIntType == "unsigned") { return strat_random<unsigned>(p_raster,numSamples,numStrata,allocation,weights,mindist,nullptr,"",0,0,filename); }
		if (minIndexIntType == "unsigned_long") { return strat_random<unsigned long>(p_raster,numSamples,numStrata,allocation,weights,mindist,nullptr,"",0,0,filename); }
		if (minIndexIntType == "unsigned_long_long") { return strat_random<unsigned long long>(p_raster,numSamples,numStrata,allocation,weights,mindist,nullptr,"",0,0,filename); }
	}
	else { //method == "Queinnec"
		if (minIndexIntType == "unsigned_short") { return strat_queinnec<unsigned short>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,nullptr,"",0,0,filename); }
		if (minIndexIntType == "unsigned") { return strat_queinnec<unsigned >(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,nullptr,"",0,0,filename); }
		if (minIndexIntType == "unsigned_long") { return strat_queinnec<unsigned long>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,nullptr,"",0,0,filename); }
		if (minIndexIntType == "unsigned_long_long") { return strat_queinnec<unsigned long long>(p_raster,numSamples,numStrata,wrow,wcol,allocation,weights,mindist,nullptr,"",0,0,filename); }
	}

}

PYBIND11_MODULE(strat, m) {
	m.def("strat_cpp", &strat_cpp);
	m.def("strat_cpp_access", &strat_cpp_access);
}
