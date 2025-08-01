/******************************************************************************
 *
 *
 * Project: sgs
 * Purpose: Integrate BalancedSampling package into sgs
 * Author: Joseph Meyer
 * Date: July, 2025
 *
 * Using the BalancedSampling package:
 * https://github.com/envisim/BalancedSampling/tree/2.1.1/src
 * which has been adapted for use directly from C++ rather than R
 * 
 * License: GPL (>=2)
 *
 ******************************************************************************/

#include <iostream>

//sgs/utils cpp code
#include "raster.h"
#include "vector.h"
#include "access.h"

//Balanced Sampling package
#include "CubeClass.h"
#include "CubeStratifiedClass.h"
#include "IndexListClass.h"
#include "KDTreeClass.h"
#include "LpmClass.h"

/**
 * This function conducts balanced sampling on the given raster image.
 *
 * To start out, the raster bands are iterated through as scanlines
 * and copied over to a vector whose data will be passed to the balanced
 * sampling function. They are copied, because the pixels being passed
 * to the BalancedSampling function should not be nan, and only the data
 * pixels are copied. 
 *
 * During this iteration two vectors keep track of specific indexes and 
 * their adjusted index after the nodata pixels are removed. 
 * These vectors are later used to calculate the original index (and thus points) 
 * of the samples. This is done instead of keeping track of x/y values
 * for every pixel, because for large images keeping track of the x/y values
 * of every pixel would be quite memory intensive.
 *
 * The probabilities for each pixel are either taken from a user input
 * (py buffer) or are set to be equal for all pixels.
 *
 * The data, as well as the probabilities and any required dimension
 * information is passed to functions from the BalancedSampling R 
 * package. The functions/classes are all in C++ and so are complied
 * and used by C++ without their R wrappers.
 *
 * Once the samples are returned, the points are calculated and returned
 * as a GDALVectorWrapper containing a vector dataset with the points.
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	std::vector<size_t> bandIndexes,
	GDALRasterWrapper *p_sraster,
	size_t stratBand,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string method,
 	py::buffer prob,
	std::string filename) 
{
	//step 1: set default paramters according to
	//https://github.com/envisim/BalancedSampling/blob/2.1.1/R/lcube.R
	//https://github.com/envisim/BalancedSampling/blob/2.1.1/R/utils.R
	size_t treeBucketSize = 50;
	double eps = 1e-12;
	int treeMethod = 2; //'kdtree2' method
	
	//step 2: determine raster info from passed raster
	size_t maxIndex = std::numeric_limits<size_t>::max();
	size_t height = static_cast<size_t>(p_raster->getHeight());
	size_t width = static_cast<size_t>(p_raster->getWidth());
	size_t bandCount = bandIndexes.size();

	//error check allocated raster size
	if (
		(maxIndex / std::max(bandCount, static_cast<size_t>(2))) < sizeof(double) ||
		(maxIndex / width) < (std::max(bandCount, static_cast<size_t>(2)) * sizeof(double)) ||
		(maxIndex / height) < (width * std::max(bandCount, static_cast<size_t>(2)) * sizeof(double))
	) {
		throw std::runtime_error("max index is too large to be processed, because the balanced sampling package requires the full raster in memory.");
	}

	if (method == "lpm2_kdtree") {
		throw std::runtime_error("this not implemented yet.");
	}

	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using band
		p_mask = maskInfo.second; //pointer to mask
	}

	//allocate spread, balanced, and potentially strata matrices, error check their allocation
	double * p_balanced = (double *) VSIMalloc3(height, width, bandCount * sizeof(double));
	double * p_spread = (double *) VSIMalloc3(height, width, 2 * sizeof(double));
	int * p_strata = (int *)((method == "lcubestratified") ? nullptr : VSIMalloc3(height, width, sizeof(int)));
	if (!p_balanced || !p_spread || (method == "lcubestratified" && !p_strata)) {
		throw std::runtime_error("unable to allocate raster band to memory");
	}	

	std::vector<double> bandNoDataValues;
	for (size_t i = 0; i < bandCount; i++) {
		GDALRasterBand *p_band = p_raster->getDataset()->GetRasterBand(bandIndexes[i] + 1);
		bandNoDataValues.push_back(p_band->GetNoDataValue());
		CPLErr err = p_band->RasterIO(
			GF_Read,
			0,
			0,
			width,
			height,
			(void *)((size_t)p_balanced + i * sizeof(double)),
			width,
			height,
			GDT_Float64,
			bandCount * sizeof(double),
			width * bandCount * sizeof(double)
		);
		if (err) {
			throw std::runtime_error("CPL error reading raster band, CPL error code: " + std::to_string(err));
		}
	}
	if (method == "lcubestratified") {
		GDALRasterBand *p_band = p_sraster->getDataset()->GetRasterBand(stratBand + 1);
		bandNoDataValues.push_back(p_band->GetNoDataValue());
		CPLErr err = p_band->RasterIO(
			GF_Read,
			0,
			0,
			width,
			height,
			p_strata,
			width,
			height,
			GDT_Int32,
			0,
			0
		);
		if (err) {
			throw std::runtime_error("CPL error reading raster band, CPL error code: " + std::to_string(err));
		}
	}

	size_t noDataCount = 0;
	size_t readIndex = 0;
	size_t writeIndex = 0;
	for (size_t y = 0; y < height; y++) {
		for (size_t x = 0; x < width; x++) {
			p_spread[writeIndex * 2] = static_cast<double>(x);
			p_spread[writeIndex * 2 + 1] = static_cast<double>(y);

			bool isNan = false;
			isNan |= (p_access && ((uint8_t *)p_mask)[readIndex] == 0);

			if (method == "lcubestratified") {
				int val = p_strata[readIndex];
				isNan |= (std::isnan(val) || val == static_cast<int>(bandNoDataValues.back()));
				if (readIndex != writeIndex) {
					p_strata[writeIndex] = val;
				}	
			}

			size_t readPixelStart = readIndex * bandCount;
			size_t writePixelStart = writeIndex * bandCount;
			for (size_t j = 0; j < bandCount; j++) {
				double val = p_balanced[readPixelStart + j];
				isNan |= (std::isnan(val) || val == bandNoDataValues[j]);
				if (readIndex != writeIndex) {
					p_balanced[writePixelStart + j] = val;
				}
			}

			readIndex++;
			writeIndex += !isNan;
			noDataCount += isNan;
		}
	}
	if (p_access) {
		free(p_accessMaskDataset);
	}

	size_t dataPixelCount = (height * width) - noDataCount;
	VSIRealloc(p_balanced, dataPixelCount * bandCount * sizeof(double));
	VSIRealloc(p_spread, dataPixelCount * 2 * sizeof(double));
	if (!p_balanced || !p_spread) {
		throw std::runtime_error("unable to re allocate rasters");
	}

	//determine prob == if prob is not already defined, ensure noData or inaccessable pixels are set to probability 0
	double *p_prob;
	std::vector<double> probVect;
	std::vector<double> xbal;
	if (prob.request().ndim != 1) {
		throw std::runtime_error("this messes up future calculation");
	}

	size_t probLength = static_cast<size_t>(prob.request().shape[0]);
	if (probLength != 0 && probLength != (height * width) - noDataCount) {
		std::cout << "**warning** length of prob list provided (" 
			  << std::to_string(probLength)
			  << ") does not match size of band (" 
			  << std::to_string(height * width) 
			  << ")." << std::endl;

		std::cout << "creating prob of correct size with equal probabilities." << std::endl;
	}

	if (probLength != height * width) {
		probVect = std::vector<double>(dataPixelCount, (double)numSamples / (double)(dataPixelCount));
		p_prob = probVect.data();
	}
	else {
		p_prob = (double *)prob.request().ptr;
	}

	//call BalancedSampling function according to desired method
	std::vector<size_t> *indexes;
	Cube * p_cube = nullptr;
	CubeStratified * p_scube = nullptr;
	if (method == "lcube") {
		p_cube = new Cube(
    			p_prob, 		//const double *
    			p_balanced, 		//double *
    			dataPixelCount,		//const size_t
    			bandCount,		//const size_t
    			eps,			//const double
    			p_spread, 		//double *
    			2,			//const size_t
    			treeBucketSize,		//const size_t
    			treeMethod		//const int
  		);
		p_cube->Run();
		indexes = &(p_cube->sample);
	}
	else if (method == "lcubestratified") {
		p_scube = new CubeStratified(
			p_strata,		//int *
			p_prob,			//const double *
			p_balanced,		//double *
			dataPixelCount,		//const size_t
			bandCount,		//const size_t
			eps,			//const double
			p_spread,		//double *
			2,			//const size_t
			treeBucketSize,		//const size_t
			treeMethod		//const int
		);
		p_scube->Run();
		indexes = &(p_scube->sample_);
	}

	//free up any allocated matrices
	CPLFree(p_balanced);
	CPLFree(p_strata);

	//step X: create GDAL Dataset to hold samples
	//TODO error check this???
	GDALAllRegister();
	GDALDataset *p_sampleDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	OGRLayer *p_sampleLayer = p_sampleDataset->CreateLayer("samples", nullptr, wkbPoint, nullptr);

	std::vector<double> xCoords;
	std::vector<double> yCoords;
	double *GT = p_raster->getGeotransform();

	//step X: turn the BalancedSampling return samples into their original indexes and calculate points to add
	for (const size_t& index : *indexes) {
		//calculate and create coordinate point
		double y = p_spread[index * 2 + 1];
		double x = p_spread[index * 2];
		double yCoord = GT[3] + x * GT[4] + y * GT[5];
		double xCoord = GT[0] + x * GT[1] + y * GT[2];
		OGRPoint newPoint = OGRPoint(xCoord, yCoord);
			
		//add point to dataset
		OGRFeature *p_feature = OGRFeature::CreateFeature(p_sampleLayer->GetLayerDefn());
		p_feature->SetGeometry(&newPoint);
		p_sampleLayer->CreateFeature(p_feature);
		OGRFeature::DestroyFeature(p_feature);

		//add to xCoords and yCoords for plotting
		xCoords.push_back(xCoord);
		yCoords.push_back(yCoord);
	}
	
	//step X: free up allocated BalancedSampling class
	free(p_spread);
	delete p_cube;
	delete p_scube;

	//step X: create new GDALVectorWrapper using dataset
	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_sampleDataset);

	//step X: write to file if filename is given
	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown while trying to write file: " << e.what() << std::endl;
		}
	}

  	return {{xCoords, yCoords}, p_sampleVectorWrapper};

	/*
	//step 3: initialize vectors to pass raster data to balancedSampling classes
	std::vector<double> xspread(height * width * bandCount);
	std::vector<int> strata;
	if (p_sraster) {
		strata = std::vector<int>(height * width);
	}	

	//step 4: allocate bands which will hold raster scanline buffers
	std::vector<double *> bands(bandCount, nullptr);
	int *p_strataBand;
	for (size_t band = 0; band < bandCount; band++) {
		bands[band] = (double *) VSIMalloc2(width, sizeof(double));
		if (!bands[band]) {
			throw std::runtime_error("unable to allocate band " + std::to_string(band + 1));
		}
	}
	if (p_sraster) {
		p_strataBand = (int *) VSIMalloc2(width, sizeof(int));
		if (!p_strataBand) {
			throw std::runtime_error("unable to allocate strat raster band");
		}
	}
	
	//step 5: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using band
		p_mask = maskInfo.second; //pointer to mask
	}

	double noDataValue = p_raster->getDataset()->GetRasterBand(1)->GetNoDataValue();
	int strataNoDataValue = -1;
	if (method == "lcubestratified") {
		strataNoDataValue = static_cast<int>(p_sraster->getDataset()->GetRasterBand(stratBand + 1)->GetNoDataValue());
		std::cout << "strataNoDataValue is " << strataNoDataValue << std::endl;
	}
	
	//step 5: initialize and set up data structures used to keep track of the original raster indexes
	//since only data pixels will be passed.
	size_t noDataCount = 0;
	std::vector<size_t> originalIndexes(width * 2);
	std::vector<size_t> adjustedIndexes(width * 2);
	size_t storedIndexesIndex = 0;
	bool prevNan = false;
	bool isNan = false;
	size_t bandIndex, originalIndex, adjustedIndex, xspreadIndex;

	//step 6: Iterate through raster scanlines and bands
	for (size_t y = 0; y < height; y++) {
		//step 6.1: read scanline from each raster band into bands vector
		for (size_t i = 0; i < bandCount; i++) {
			CPLErr err = p_raster->getDataset()->GetRasterBand(bandIndexes[i] + 1)->ReadRaster(
				bands[i],			//pData
			       	width,				//nArrayEltCount
				0,				//dfXOff
				y,				//dfYOff
				width,				//dfXSize
				1,				//dfYSize
				width,				//nBufXSize
				1,				//nBufYSize
				GRIORA_NearestNeighbour,	//eResampleAlg (default)
				nullptr,			//pfnProgress
				nullptr				//pProgressData
			);
			if (err) {
				std::string errorStr = "ReadRaster failed with CPLError code: " + std::to_string(err);
				throw std::runtime_error(errorStr);
			}
		}
		if (p_sraster) {
			CPLErr err = p_sraster->getDataset()->GetRasterBand(stratBand + 1)->ReadRaster(
				p_strataBand,			//pData
				width,				//nArrayEltCount
				0,				//dfXOff
				y,				//dfYOff
				width,				//dfXSize
				1,				//dfYSize
				width,				//nBufXSize
				1,				//nBufYSize
				GRIORA_NearestNeighbour,	//eResampleAlg (default)
				nullptr,			//pfnProgress
				nullptr				//pProgressData
			);
			std::string errorStr = "ReadRaster failed with CPLError code: " + std::to_string(err);
			throw std::runtime_error(errorStr);

		}

		size_t startIndex = y * width;
		
		//step 6.2: special case for first index, to set up originalIndexes and adjustedIndexes vectors as required
		if (y == 0) {
			for (size_t i = 0;i < bands.size(); i++) {
				double val = bands[i][0];
				xspread[i] = val;
				isNan |= (std::isnan(val) || val == noDataValue);
			}
			if (p_sraster) {
				int val = p_strataBand[0];
				strata[0] = val;
				isNan |= val == strataNoDataValue;
			}
			isNan |= (p_access && ((uint8_t *)p_mask)[0] == 0);
			prevNan = isNan;
			noDataCount += isNan;

			//update the adjusted and original indexes, they will be overwritten if isNan is true
			originalIndexes[storedIndexesIndex] = 0;
			adjustedIndexes[storedIndexesIndex] = 0;
			storedIndexesIndex += (!isNan);
		}

		//step 6.3: iterate through raster bands copying to vector and keeping track of nodata
		for (size_t x = 0 + (y == 0); x < width; x++) {
			isNan = false;
			bandIndex = x;
			xspreadIndex = (startIndex + x - noDataCount) * bandCount;
			originalIndex = startIndex + x;
			adjustedIndex = startIndex + x - noDataCount;
			for (size_t i = 0; i < bandCount; i++) {
				double val = bands[i][bandIndex];
				xspread[xspreadIndex + i] = val;
				isNan |= (std::isnan(val) || val == noDataValue);
			}	
			if (p_sraster) { 
				int val = p_strataBand[bandIndex];
				strata[adjustedIndex] = val;
				isNan |= val == strataNoDataValue;
			}
			isNan |= (p_access && ((uint8_t *)p_mask)[originalIndex] == 0);

			//update the adjusted and original indexes, they will be overwritten if isNan is true	
			originalIndexes[storedIndexesIndex] = originalIndex;
			adjustedIndexes[storedIndexesIndex] = adjustedIndex;
			storedIndexesIndex += (!isNan && prevNan);
			
			noDataCount += isNan;
			prevNan = isNan;
		}

		//resize vectors as required
		xspread.resize(height * width * bandCount - noDataCount * bandCount);
		if (originalIndexes.size() < storedIndexesIndex + (width / 2)) {
			originalIndexes.resize(originalIndexes.size() + width * 4);
			adjustedIndexes.resize(originalIndexes.size() + width * 4);
		}
	}
	if (p_access) {
		free(p_accessMaskDataset);
	}

	size_t noNanBandSize = height * width - noDataCount;

	//step 7: free no longer used band data and adjust size of vectors
	for (size_t band = 0; band < bands.size(); band++) {
		CPLFree(bands[band]);
	}
	xspread.resize(noNanBandSize * bandCount);
	xspread.shrink_to_fit();
	originalIndexes.resize(storedIndexesIndex + 1);
	originalIndexes.shrink_to_fit();
	adjustedIndexes.resize(storedIndexesIndex + 1);
	originalIndexes.shrink_to_fit();
	if (p_sraster) {
		strata.resize(noNanBandSize);
		strata.shrink_to_fit();
	}

	//step 8: determine probability list (vector)
	double *p_prob;
	std::vector<double> probVect, xbal;

	//TODO remove
	if (prob.request().ndim != 1) {
		throw std::runtime_error("this messes up future calculation");
	}

	size_t probLength = static_cast<size_t>(prob.request().shape[0]);
	if (probLength != 0 && probLength != noNanBandSize) {
		std::cout << "**warning** length of prob list provided (" 
			  << std::to_string(probLength)
			  << ") does not match size of band without nan values (" 
			  << std::to_string(noNanBandSize) 
			  << ")." << std::endl;

		std::cout << "creating prob of correct size with equal probabilities." << std::endl;
	}

	if (probLength != noNanBandSize) {
		probVect = std::vector<double>(noNanBandSize, (double)numSamples / (double)noNanBandSize);
		p_prob = probVect.data();
	}
	else {
		p_prob = (double *)prob.request().ptr;
	}

	if (method != "lpm_kdtree") {
		//lpm_kdtree doesn't need xbal, but the others do
		xbal.resize(noNanBandSize);
		std::memcpy(
			reinterpret_cast<void *>(xbal.data()),	//dst
			reinterpret_cast<void *>(p_prob), 	//src
			noNanBandSize * sizeof(double)		//num bytes
		);
	}
	
	//step 9: call BalancedSampling function according to desired method
	std::vector<size_t> *samples;
	Cube * p_cube = nullptr;
	CubeStratified * p_scube = nullptr;
	Lpm * p_lpm = nullptr;
	if (method == "lcube") {
		p_cube = new Cube(
    			p_prob, 		//const double*
    			xbal.data(), 		//double *
    			noNanBandSize,		//const size_t
    			1,			//const size_t
    			eps,			//const double
    			xspread.data(), 	//double *
    			bandCount,		//const size_t
    			treeBucketSize,		//const size_t
    			treeMethod		//const int
  		);
		p_cube->Run();
		samples = &(p_cube->sample);
	}
	else if (method == "lcubestratified") {
		p_scube = new CubeStratified(
			strata.data(),		//int *
			p_prob,			//const double *
			xbal.data(),		//double *
			noNanBandSize,		//const size_t
			1,			//const size_t
			eps,			//const double
			xspread.data(),		//double *
			bandCount,		//const size_t
			treeBucketSize,		//const size_t
			treeMethod		//const int
		);
		p_scube->Run();
		samples = &(p_scube->sample_);
	}
	else {// method == "lpm2_kdtree"
		p_lpm = new Lpm(
			LpmMethod::lpm2,	//const LpmMethod
			p_prob,			//const double *
			xspread.data(),		//double *
			noNanBandSize,		//const size_t
			bandCount,		//const size_t
			eps,			//const double
			treeBucketSize,		//const size_t
			treeMethod		//const int
		);
		p_lpm->Run();
		samples = &(p_lpm->sample);
	}
 
	//step 10: allow memory to be freed if helpful
	xspread.clear();
	xspread.shrink_to_fit();
	strata.clear();
	strata.shrink_to_fit();
	probVect.clear();
	probVect.shrink_to_fit();
	xbal.clear();
	xbal.shrink_to_fit();

	//step 11: create GDAL Dataset to hold samples
	//TODO error check this???
	GDALAllRegister();
	GDALDataset *p_sampleDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	OGRLayer *p_sampleLayer = p_sampleDataset->CreateLayer("samples", nullptr, wkbPoint, nullptr);

	std::vector<double> xCoords;
	std::vector<double> yCoords;
	double *GT = p_raster->getGeotransform();

	//step 12: turn the BalancedSample return samples into their original indexes and calculate points to add
	for (const size_t& sample : *samples) {
		//get original index from adjusted sample index using the adjustedIndexes and originalIndexes vectors
		auto boundIt = std::upper_bound(adjustedIndexes.begin(), adjustedIndexes.end(), sample);
		size_t boundIndex = (boundIt == adjustedIndexes.end()) ? adjustedIndexes.back() : std::distance(adjustedIndexes.begin(), boundIt) - 1;
		size_t index = sample + originalIndexes[boundIndex] - adjustedIndexes[boundIndex];

		//calculate and create coordinate point
		double yIndex = index / p_raster->getWidth();
		double xIndex = index - (yIndex * p_raster->getWidth());
		double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
		double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
		OGRPoint newPoint = OGRPoint(xCoord, yCoord);
			
		//add point to dataset
		OGRFeature *p_feature = OGRFeature::CreateFeature(p_sampleLayer->GetLayerDefn());
		p_feature->SetGeometry(&newPoint);
		p_sampleLayer->CreateFeature(p_feature);
		OGRFeature::DestroyFeature(p_feature);

		//add to xCoords and yCoords for plotting
		xCoords.push_back(xCoord);
		yCoords.push_back(yCoord);
	}

	//step 13: free up pointers which may have been allocated
	delete p_cube;
	delete p_scube;
	delete p_lpm;

	//step 14: create new GDALVectorWrapper using dataset
	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_sampleDataset);

	//step 15: write to file if filename is given
	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown while trying to write file: " << e.what() << std::endl;
		}
	}

  	return {{xCoords, yCoords}, p_sampleVectorWrapper};
	*/
}

/**
 * balanced sampling function which is called when the user does
 * not include either an access vector or a strat raster.
 *
 * calls balanced() where all access-related and sraster-related
 * paramters are null/0/"".
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced_cpp(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	std::vector<size_t> bandIndexes,
	std::string method,
	py::buffer prob,
	std::string filename)
{
	return balanced(
		p_raster,
		numSamples,
		bandIndexes,
		nullptr,
		0,
		nullptr,
		"",
		0,
		0,
		method,
		prob,
		filename
	);
}

/**
 * balanced sampling function when the user does not include a strat
 * raster.
 *
 * calls balanced() where all sraster-related parameters are
 * null/0/"".
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced_access_cpp(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	std::vector<size_t> bandIndexes,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string method,
	py::buffer prob,
	std::string filename)
{
	return balanced(
		p_raster,
		numSamples,
		bandIndexes,
		nullptr,
		0,
		p_access,
		layerName,
		buffInner,
		buffOuter,
		method,
		prob,
		filename
	);
}

/**
 * balanced sampling function when the user does not include an
 * access vector.
 *
 * calls balanced() where all access-related parameters are null/0/"".
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced_strata_cpp(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	std::vector<size_t> bandIndexes,
	GDALRasterWrapper *p_sraster,
	size_t stratBand,
	std::string method,
	py::buffer prob,
	std::string filename)
{
	return balanced(
		p_raster,
		numSamples,
		bandIndexes,
		p_sraster,
		stratBand,
		nullptr,
		"",
		0,
		0,
		method,
		prob,
		filename
	);
}


PYBIND11_MODULE(balanced, m) {
	m.def("balanced_cpp", &balanced_cpp);
	m.def("balanced_access_cpp", &balanced_access_cpp);
	m.def("balanced_strata_cpp", &balanced_strata_cpp);
	m.def("balanced_access_strata_cpp", &balanced);
}
