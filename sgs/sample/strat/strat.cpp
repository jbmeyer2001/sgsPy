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
typedef <typename U>
std::vector<size_t>
calculateAllocation(
	size_t numSamples,
	std::string allocation, 
	std::vector<size_t>& strataSizes, 
	std::vector<double>& weights,
	size_t numPixels
{
	std::vector<U> retval;

	if (allocation == "prop") {
		//allocate the samples per stratum according to stratum size
		size_t remainer = numSamples;
		for (size_t i = 0; i < strataSizes.size(); i++) {
			size_t count = numSamples / (numPixels / strataSizes[i]);
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
		size_t strataSampleCount = numStrata / stratSizes.size();
		
		for (size_t i = 0; i < strataSizes.size(); i++) {
			if (strataSizes[i] < strataSampleCount) {
				retval.push_back(strataSizes[i]);
				std::cout << "warning: strata " << i << " does not have enough pixels for the full " << strataSampleCount << " samples it should recieve. There will be less than " << numSamples << " final samples." << std::endl;
			}
			else {
				retval.push_back(strataSampleCount);
			}

		}
	}
	else if (allocation == "manual") { //TODO what if a strata doesn't have enough pixels???
		if (strataSizes.size() != weights.size()) {
			throw std::runtime_error("size of weights is not equal to the number of strata.");
		}

		//allocate samples accordign to weights.
		remainder = numSamples;
		for (size_t i = 0; i < strataSizes.size(); i++) {
			size_t count = (size_t)((double)strataSizes[i] * weights[i]);
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
			if (retval[i] > strataSizes[i]) {
				std::cout << "warning: strata " << i << " does not have enough pixels for the full " << retval[i] << " samples it should recieve. There will be less than " << numSamples << " final samples." << std::endl;
				retval[i] = strataSizes[i];
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
	return isnan;
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
	float *p_strat = p_raster->getBand(0);

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
	double noDataValue = p_dataset->GetDataset()->GetRasterBand(1)->GetNoDataValue();
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
		strataSizes,
		weights,
		numDataPixels
	);
	U maxSampleCount = *std::max_element(stratumCounts.begin(), stratumCounts.end());
	
	//Step 6: generate random number generator using mt19937	
	std::mt19937 gen(time(nullptr));
	std::vector<std::function<U(void)>> rngs;
	for (size_t i = 0; i < stratumCounts.size(); i++) {
		rngs.push_back(std::bind(
			std::ref(std::uniform_int_distribution<U>(0, stratumCounts[i] - 1)),
			std::ref(gen())
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
			while(backupSamplePixelsSize() < startSize + (stratumCounts[i] * 2)) {
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
	std::Vector<std::string> wktPoints;

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
std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
strat_queinnec(
	GDALRasterWrapper,
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
	//add later	
	return {{}, ""};
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
