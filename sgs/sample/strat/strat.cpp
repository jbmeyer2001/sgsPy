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

/*
 *
 */
template <typename U>
inline bool checkNan(U index, float *val, uint8_t *p_mask, float *p_strata, double noDataValue) {
	bool isNan = p_mask && p_mask[index] == 0;
	if (!isNan) {
		*val = p_strata[index];
		isNan = std::isnan(*val) || static_cast<double>(*val) == noDataValue;
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

	//step 7: determine pixel values to include as samples.
	
	//make a set of pixels to include as samples
	std::vector<std::unordered_set<U>> sampleIndexes;
	std::vector<typename std::unordered_set<U>::iterator> sampleIterators;
	std::vector<U> samplesAdded;
	std::vector<U> strataNum;

	std::mt19937::result_type seed = time(nullptr);
	for (size_t i = 0; i < stratumCounts.size(); i++) {	
		sampleIndexes.push_back({});
		samplesAdded.push_back(0);
		strataNum.push_back(i);

		auto rng = std::bind(
			std::uniform_int_distribution<U>(0, stratumIndexes[i].size() - 1),
			std::mt19937(seed)
		);

		U stratumSamples = std::min((mindist == 0) ? (size_t)stratumCounts[i] : (size_t)stratumCounts[i] * 3, stratumIndexes[i].size());

		if (stratumSamples > stratumIndexes[i].size() / 2) {
			std::unordered_set<U> dontSamplePixels;
			while (dontSamplePixels.size() < stratumIndexes[i].size() - stratumSamples) {
				dontSamplePixels.insert(rng());
			}
			for (size_t j = 0; j < stratumIndexes[i].size(); j++) {
				if (dontSamplePixels.find(j) == dontSamplePixels.end()) {
					sampleIndexes[i].insert(stratumIndexes[i][j]);
				}
			}
		}
		else {
			size_t samplePixelsStartSize = sampleIndexes.size();
			while (sampleIndexes.size() < samplePixelsStartSize + stratumSamples) {
				sampleIndexes[i].insert(stratumIndexes[i][rng()]);
			}
		}

		sampleIterators.push_back(sampleIndexes[i].begin());
	}

	//step 8: generate coordinate points for each sample index, and only add if they're outside of mindist
	std::vector<double> xCoords;
	std::vector<double> yCoords;
	std::vector<OGRPoint> points;
	std::vector<std::string> wktPoints;

	U sIndex = 0;
	U completedStratum = 0;
	double *GT = p_raster->getGeotransform();
	while (completedStratum < stratumCounts.size()) {
		//determine strata from i
		if (sIndex >= strataNum.size()) {
			sIndex = 0;
		}
		U strata = strataNum[sIndex];

		if (sampleIterators[strata] == sampleIndexes[strata].end() || stratumCounts[strata] == samplesAdded[strata]) {
			completedStratum++;
			strataNum.erase(strataNum.begin() + sIndex);
		}
		else {
			U index = *sampleIterators[strata];
			sampleIterators[strata] = std::next(sampleIterators[strata]);

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
					sIndex++;
					continue;
				}
			}

			samplesAdded[strata]++;
			points.push_back(newPoint);
			wktPoints.push_back(newPoint.exportToWkt());
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
			sIndex++;
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
	std::string filename)
{
	//step 1: get raster band
	if (p_raster->getRasterType() != GDT_Float32) {
		throw std::runtime_error("raster band must be of type float32.");
	}
	float *p_strata = (float *)p_raster->getRasterBand(0);

	//step 2: create stratum index storing vectors
	std::vector<std::vector<U>> queinnecStratumIndexes;
	std::vector<std::vector<U>> randomStratumIndexes;
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
	U width = p_raster->getWidth();
	U height = p_raster->getHeight();
	U horizontalPad = wcol / 2;
	U verticalPad = wcol / 2;
	U fwHeight = wrow;
	U fwWidth = width - wcol + 1;
	std::vector<bool> focalWindowMatrix(fwWidth * fwHeight, true);
	std::vector<bool> prevVertSame(width, true);
	bool prevHoriSame = true;
	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
	U noDataPixelCount = 0;
	bool nextVertSame;
	bool nextHoriSame;

	U y = 0;
	int64_t fwy = -wrow + 1;
	U x;
	int64_t fwx;
	bool addSelf;
	bool addfw;

	while (y < height) {
		//reset no longer used section of focal window matrix
		for (int64_t fwxi = 0; fwxi < fwWidth; fwxi++) {
			focalWindowMatrix[((fwy - 1) % wrow) * fwWidth + fwxi] = true;
		}

		addSelf = (fwy + verticalPad < 0) || (y - verticalPad > height - wrow + 1);
		addfw = fwy >= 0;

		int64_t fwyStart = std::max(fwy, static_cast<int64_t>(0));
		int64_t fwyEnd = std::min(y, static_cast<U>(height - wrow + 1));
		
		x = 0;
		fwx = -wcol + 1;

		while (x < width) {
			addSelf |= (fwx + horizontalPad < 0) || (x - verticalPad > fwWidth);
			addfw &= fwx >= 0;

			int64_t fwxStart = std::max(fwx, static_cast<int64_t>(0));
			int64_t fwxEnd = std::min(x, fwWidth);
		
			U index = y * width + x;
			float val;
			bool isNan = checkNan<U>(index, &val, p_mask, p_strata, noDataValue);
			noDataPixelCount += (U)isNan;

			nextVertSame = !isNan && ((y == height - 1) || val == p_strata[index + width]);
		       	nextHoriSame = !isNan && ((x == width - 1) || val == p_strata[index + 1]);	

			if (addSelf && !isNan) {
				randomStratumIndexes[(size_t)val].push_back(index);
			}

			if (addfw) {
				int64_t fwIndex = (fwy % wrow) * fwWidth + fwx;
				index = (fwy + verticalPad) * width + fwx + verticalPad;	
				bool isNan = checkNan<U>(index, &val, p_mask, p_strata, noDataValue);

				if (!isNan && focalWindowMatrix[fwIndex]) {
					queinnecStratumIndexes[(size_t)val].push_back(index);
				}
				else if (!isNan) {
					randomStratumIndexes[(size_t)val].push_back(index);
				}

			}

			int64_t fwi;
			//set the portion of the focal window matrix to false which is impacted
			//only when the next horizontal pixel value is different than the 
			//current one.
			if (!nextHoriSame) {
				for (int64_t fwyi = fwyStart + 1; fwyi <= fwyEnd - 1; fwyi++) {
					fwi = (fwyi % wrow) * fwWidth + fwxEnd;
					focalWindowMatrix[fwi] = false;
				}
			}

			//set the portion of the focal window matrix to false which is impacted
			//only when the next vertical pixel value is different than the current
			//one.
			if (!nextVertSame) {
				for (int64_t fwxi = fwxStart + 1; fwxi <= fwxEnd - 1; fwxi++) {
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
				for (int64_t fwxi = fwxStart + 1; fwxi <= fwxEnd - 1; fwxi++) {
					fwi = (fwyStart % wrow) * fwWidth + fwxi;
					focalWindowMatrix[fwi] = false;
				}
			}

			//set the portion of the focal window matrix to false which must be
			//changed when the next vertical pixel is different but the previous
			//vertical pixel is the same as teh current one. 
			if (!nextVertSame && prevVertSame[x]) {
				for (int64_t fwyi = fwyStart + 1; fwyi <= fwyEnd - 1; fwyi++) {
					fwi = (fwyi % wrow) * fwWidth + fwxStart;
					focalWindowMatrix[fwi] = false;
				}
			}

			//set the portion of the focal window matrix to false which must be
			//changed when either the next vertical or horizontlal pixel is different,
			//but both the previous vertical and previous horizontal pixels were the same
			//as the current one.
			if ((!nextHoriSame || !nextVertSame) && prevHoriSame && prevVertSame[x]) {
				for (int64_t fwyi = fwyStart + 1; fwyi <= fwyEnd - 1; fwyi++) {
					for (int64_t fwxi = fwxStart + 1; fwxi <= fwxEnd - 1; fwxi++) {
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

	//step 8: calculate allocation of samples depending on stratum sizes 
	std::vector<U> strataSizes;
	for (size_t i = 0; i < queinnecStratumIndexes.size(); i++) {
		strataSizes.push_back(randomStratumIndexes[i].size() + queinnecStratumIndexes[i].size());
	}
	U numDataPixels = p_raster->getWidth() * p_raster->getHeight() - noDataPixelCount;

	std::vector<U> stratumCounts = calculateAllocation<U>(
		numSamples,
		allocation,
		&strataSizes,
		&weights,
		numDataPixels
	);

	//step 7: determine queinnec indexes to try including as samples
	std::vector<std::unordered_set<U>> sampleIndexes;
	std::vector<typename std::unordered_set<U>::iterator> sampleIterators;
	std::vector<U> samplesAdded;
	std::vector<U> strataNum;

	std::mt19937::result_type seed = time(nullptr);
	for (size_t i = 0; i < stratumCounts.size(); i++) {	
		sampleIndexes.push_back({});
		samplesAdded.push_back(0);
		strataNum.push_back(i);

		auto rng = std::bind(
			std::uniform_int_distribution<U>(0, queinnecStratumIndexes[i].size() - 1),
			std::mt19937(seed)
		);

		U stratumSamples = std::min((mindist == 0) ? (size_t)stratumCounts[i] : (size_t)stratumCounts[i] * 3, queinnecStratumIndexes[i].size());

		if (stratumSamples > queinnecStratumIndexes[i].size() / 2) {
			std::unordered_set<U> dontSamplePixels;
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
			size_t samplePixelsStartSize = sampleIndexes.size();
			while (sampleIndexes.size() < samplePixelsStartSize + stratumSamples) {
				sampleIndexes[i].insert(queinnecStratumIndexes[i][rng()]);
			}
		}

		sampleIterators.push_back(sampleIndexes[i].begin());
	}

	//step 8: generate coordinate points for each queinnec sample, and only add if they're outside of mindist
	std::vector<double> xCoords;
	std::vector<double> yCoords;
	std::vector<OGRPoint> points;
	std::vector<std::string> wktPoints;

	U sIndex = 0;
	U completedStratum = 0;
	double *GT = p_raster->getGeotransform();
	while (completedStratum < stratumCounts.size()) {
		//determine strata from i
		if (sIndex >= strataNum.size()) {
			sIndex = 0;
		}
		U strata = strataNum[sIndex];

		if (sampleIterators[strata] == sampleIndexes[strata].end() || stratumCounts[strata] == samplesAdded[strata]) {
			completedStratum++;
			strataNum.erase(strataNum.begin() + sIndex);
		}
		else {
			U index = *sampleIterators[strata];
			sampleIterators[strata] = std::next(sampleIterators[strata]);

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
					sIndex++;
					continue;
				}
			}

			samplesAdded[strata]++;
			points.push_back(newPoint);
			wktPoints.push_back(newPoint.exportToWkt());
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
			sIndex++;
		}
	}

	//step 9: determine random indexes to try including as samples
	for (size_t i = 0; i < stratumCounts.size(); i++) {	
		stratumCounts[i] -= samplesAdded[i];
		if (stratumCounts[i] > 0) {
			completedStratum--;
			sampleIndexes[i] = {};
			samplesAdded[i] = 0;
			strataNum.push_back(i);

			auto rng = std::bind(
				std::uniform_int_distribution<U>(0, randomStratumIndexes[i].size() - 1),
				std::mt19937(seed)
			);

			U stratumSamples = std::min((mindist == 0) ? (size_t)stratumCounts[i] : (size_t)stratumCounts[i] * 3, randomStratumIndexes[i].size());

			if (stratumSamples > randomStratumIndexes[i].size() / 2) {
				std::unordered_set<U> dontSamplePixels;
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
				size_t samplePixelsStartSize = sampleIndexes.size();
				while (sampleIndexes.size() < samplePixelsStartSize + stratumSamples) {
					sampleIndexes[i].insert(randomStratumIndexes[i][rng()]);
				}
			}

			sampleIterators[i] = sampleIndexes[i].begin();
		}
	}

	//step 10: try adding random samples
	sIndex = 0;
	while (completedStratum < stratumCounts.size()) {
		//determine strata from i
		if (sIndex >= strataNum.size()) {
			sIndex = 0;
		}
		U strata = strataNum[sIndex];

		if (sampleIterators[strata] == sampleIndexes[strata].end() || stratumCounts[strata] == samplesAdded[strata]) {
			completedStratum++;
			strataNum.erase(strataNum.begin() + sIndex);
		}
		else {
			U index = *sampleIterators[strata];
			sampleIterators[strata] = std::next(sampleIterators[strata]);

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
					sIndex++;
					continue;
				}
			}

			samplesAdded[strata]++;
			points.push_back(newPoint);
			wktPoints.push_back(newPoint.exportToWkt());
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
			sIndex++;
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
