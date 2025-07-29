/******************************************************************************
 *
 * Project: sgs
 * Purpose: Integrate BalancedSampling package into sgs
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 * Adapted from code authored by Wilmer Prentius in the following files:
 * https://github.com/envisim/BalancedSampling/blob/2.1.1/src/cube.cc
 * https://github.com/envisim/BalancedSampling/blob/2.1.1/src/cube_stratified.cc
 * https://github.com/envisim/BalancedSampling/blob/2.1.1/src/hlpm2.cc
 * License: GPL (>=2)
 *
 ******************************************************************************/

#include <iostream>

//sgs/utils cpp code
#include "raster.h"
#include "vector.h"

//Balanced Sampling package
#include "CubeClass.h"
#include "CubeStratifiedClass.h"
#include "IndexListClass.h"
#include "KDTreeClass.h"
#include "LpmClass.h"

/**
 *
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced(
	GDALRasterWrapper *p_raster,
	std::vector<size_t> bandIndexes;
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
	//set default parameters according to 
	//https://github.com/envisim/BalancedSampling/blob/2.1.1/R/lcube.R
	//https://github.com/envisim/BalancedSampling/blob/2.1.1/R/utils.R
	size_t treeBucketSize = 50;
	double eps = 1e-12;
	int treeMethod = 2; //'kdtree2' method
	
	size_t maxIndex = std::numeric_limits<size_t>::max();
	size_t height = static_cast<size_t>(p_raster->getHeight());
	size_t width = static_cast<size_t>(p_raster->getWidth());
	size_t bandCount = bandIndexes.size();

	//error checking for potential overflow problems
	//since the balanced sampling package requires that all raster data be allocated at once
	if (
		(maxIndex / bandCount) < sizeof(double) ||
		(maxIndex / width) < (bandCount * sizeof(double) ||
		(maxIndex / height) < (width * bandCount * sizeof(double)
	) {
		throw runtime_error("max index is too large to be processed by the current operating system, "
				+ "because the balanced sampling package requires the full raster to be in memory at once.");
	}

	//create vectors which will hold data pixels of bands
	std::vector<double> xspread(height * width * bandCount);
	std::vector<int> strata;
	if (p_sraster) {
		strata = std::vector<int>(height * width);
	}	

	//create band storing data structure and allocate memory
	std::vector<double *> bands(bandCount, nullptr);
	int *p_strataBand;
	for (band = 0; band < bandCount; band++) {
		bands[band] = (double *) VSIMalloc2(width, sizeof(double));
		if (!bands[band]) {
			throw std::runtime_error("unable to allocate band " + std::to_string(band + 1));
		}
	}
	if (p_sraster) {
		p_strataBand = (int *) VSIMalloc(width, sizeof(int));
	}
	
	//get access mask if defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using band
		p_mask = maskInfo.second; //pointer to mask
	}

	double noDataValue = p_raster->getDataset->GetRasterBand(1)->getNoDataValue();
	int strataNoDataValue = static_cast<int>(p_sraster->getDataset->GetRasterBand(stratBand + 1)->getNoDataValue());

	std::cout << "strataNoDataValue is " << strataNoDataValue << std::endl;
	
	//noDataCount, originalIndexes vector, and adjustedIndexes vector work to 
	//map the values of a raster band with all the noData pixels removed back
	//to it's original index.
	size_t noDataCount = 0;
	std::vector<size_t> originalIndexes(width * 2);
	std::vector<size_t> adjustedIndexes(width * 2);
	size_t storedIndexesIndex = 0;
	bool prevNan = false;
	bool isNan = false;
	size_t bandIndex, originalIndex, adjustedIndex, xspreadIndex;

	for (size_t y = 0; y < height; y++) {
		//read scanline from each raster band into bands vector
		for (size_t i = 0; i < bandCount; i++) {
			CPLErr err = p_raster->getDataset->GetRasterBand(bandIndexes[i] + 1)->ReadRaster(
				bands[band],	//pData
			       	width,		//nArrayEltCount
				0,		//dfXOff
				y,		//dfYOff
				width,		//dfXSize
				1,		//dfYSize
				width,		//nBufXSize
				1,		//nBufYSize
				nullptr,	//eResampleAlg
				nullptr,	//pfnProgress
				nullptr		//pProgressData
			);
			if (err) {
				std::cout << "ReadRaster failed with CPLError code: " << err << std::endl;
			}
		}
		if (p_sraster) {
			CPLErr err = p_sraster->getDataset->GetRasterBand(stratBand + 1)->ReadRaster(
				p_strataBand,	//pData
				width,		//nArrayEltCount
				0,		//dfXOff
				y,		//dfYOff
				width,		//dfXSize
				1,		//dfYSize
				width,		//nBufXSize
				1,		//nBufYSize
				nullptr,	//eResampleAlg
				nullptr,	//pfnProgress
				nullptr		//pProgressData
			)	
		}


		size_t startIndex = y * width;
		
		//special case for first index, to set up originalIndexes and adjustedIndexes vectors as required
		if (y == 0) {
			for (size_t band = 0; band < bands.size(); band++) {
				double val = bands[band][0];
				xspread[band] = val;
				isNan |= (std::isnan(val) || val == noDataValue);
			}
			if (p_sraster) {
				int val = p_strataBand[0];
				strata[0] = val;
				isNan |= val == strataNoDataValue;
			}
			isNan |= (p_access && ((uint8_t *)p_mask)[0] == 0)
			prevNan = isNan;
			noDataCount += isNan;

			//update the adjusted and original indexes, they will be overwritten if isNan is true
			originalIndexes[storedIndexesIndex] = 0;
			adjustedIndexes[storedIndexesIndex] = 0;
			storedIndexesIndex += (!isNan);
		}

		//iterate through indexes
		//write the values over to xspread
		//if any band as nodata all values written will be overwritten at next iteration 
		for (x = 0 + (y == 0); x < width; x++) {
			isNan = false;
			bandIndex = x;
			xspreadIndex = (startIndex + x - noDataCount) * bandCount;
			originalIndex = startIndex + x;
			adjustedIndex = startIndex + x - noDataCount;
			for (size_t band = 0; band < bandCount; band++) {
				double val = bands[band][bandIndex];
				xspread[xspreadIndex + band] = val;
				isNan |= (std::isnan(val) || val == noDataValue);
			}	
			if (p_sraster) { 
				int val = p_strataBand[bandIndex];
				strata[adjustedIndex] = val;
				isNan |= val == strataNoDataValue;
			}
			isNan = (p_access && ((uint8_t *)p_mask)[originalIndex] == 0);
		
			//update the adjusted and original indexes, they will be overwritten if isNan is true	
			originalIndexes[storedIndexesIndex] = originalIndex;
			adjustedIndexes[storedIndexesIndex] = adjustedIndex;
			storedIndexesIndex += (!isNan && prevNan);
			
			noDataCount += isNan;
			prevNan = isNan;
		}

		xspread.resize(height * width * bandCount - noDataCout * bandCount);
		
		if (originalIndexes.size() < storedIndexesIndex + (width / 2)) {
			originalIndexes.resize(originalIndexes.size() + width * 4);
			adjustedIndexes.resize(originalIndexes.size() + with * 4);
		}
	}

	size_t noNanBandSize = height * width - noDataCount;

	//free no longer used band data and adjust size of vectors
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

	//determine probability list (vector)
	double *p_prob, p_xbal;
	std::vector<double> prob, xbal;

	//TODO remove
	if (prob.request().ndim != 1) {
		throw std::runtime_error("this messes up future calculation");
	}

	ssize_t probLength = prob.request().shape[0];
	if (probLength != 0 && probLength != noNanBandSize) {
		std::cout << "**warning** length of prob list provided (" 
			  << std::to_string(probLength) +
			  << ") does not match size of band without nan values (" 
			  << std::to_string(noNanBandSize) 
			  << ")." << std::endl;

		std::cout << "creating prob of correct size with equal probabilities." << std::endl;
	}

	if (probLength != noNanBandSize) {
		prob = std::vector<double>(noDataBandSize, (double)numSamples / (double)noNanBandSize);
		p_prob = prob.data();
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
		p_xbal = xbal.data();
	}

	std::vector<size_t> *samples;
	std::unique_ptr<Cube *> p_cube;
	std::unique_ptr<CubeStratified *> p_scube;
	std::unique_ptr<Lpm *> p_lpm;
	if (method == "lcube") {
		p_cube = std::make_unique<Cube *>(
    			p_prob, 		//const double*
    			p_xbal, 		//double *
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
		p_scube = std::make_unique<CubeStratified *>(
			strata.data(),		//int *
			p_prob,			//const double *
			p_xbal,			//double *
			noNanBandSize,		//const size_t
			1,			//const size_t
			eps,			//const double
			xspread.data(),		//double *
			bandCount,		//const size_t
			treeBucketSize,		//const size_t
			treeMethod		//const int
		);
		p_scube.Run();
		samples = &(p_scube->sample_);
	}
	else {// method == "lpm2_kdtree"
		p_lpm = std::make_unique<Lpm *>(
			lpmMethod::lpm2,	//const LpmMethod
			p_prob,			//const double *
			xspread.data(),		//double *
			noNanBandSize,		//const size_t
			1,			//const size_t
			eps,			//const double
			treeBucketSize,		//const size_t
			treeMethod		//const int
		);
		p_lpm->Run();
		samples = &(p_lpm->sample);
	}
  
	//don't keep vectors full if we're not using them
	//as it may allow memory to be released
	xspread.clear();
	xspread.shrink_to_fit();
	strata.clear();
	strata.shrink_to_fit();
	prob.clear();
	prob.shrink_to_fit();
	xbal.clear();
	xbal.shrink_to_fit();
	originalIndexes.clear();
	originalIndexes.shrink_to_fit();
	adjustedIndexes.clear();
	adjustedIndexes.shrink_to_fit();

	//TODO error check this???
	GDALAllRegister();
	GDALDataset *p_sampleDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	OGRLayer *p_sampleLayer = p_sampleDataset->CreateLayer("samples", nullptr, wkbPoint, nullptr);

	std::vector<double> xCoords;
	std::vector<double> yCoords;
	double *GT = p_raster->getGeotransform();

	for (const size_t& sample : *samples) {
		//get original index from adjusted sample index using the adjustedIndexes and originalIndexes vectors
		auto boundIt = std::upper_bound(adjustedIndexes.begin(), adjustedIndexes.end(), sample);
		aize_t boundIndex = (boundIt == adjustedIndexes.end()) ? adjustedIndexes.back() : std::distance(adjustedIndexes.begin(), boundIt) - 1;
		size_t index = sample + originalIndexes[boundIndex] - adjustedIndexes[boundIndex];

		//calculate and create coordinate point
		double yIndex = index / p_raster->getWidth();
		double xIndex = index - (yIndex * p_raster->getWidth());
		double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
		double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
		OGRPoint newPoint = OGRPoint(xCoord, yCoord);
			
		//add point to dataset
		OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
		p_feature->SetGeometry(&newPoint);
		p_layer->CreateFeature(p_feature);
		OGRFeature::DestroyFeature(p_feature);

		//add to xCoords and yCoords for plotting
		xCoords.push_back(xCoord);
		yCoords.push_back(yCoord);
	}

	//free up unique_ptr which may have been allocated
	p_cube.reset();
	p_scube.reset();
	p_lpm.reset();

	//create new GDALVectorWrapper using dataset
	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_sampleDataset);

	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch () {
			std::cout << "Exception thrown while trying to write file: " << e.what() << std::endl;
		}
	}
  	return {{xCoords, yCoords}, p_sampleVectorWrapper};
}

/**
 *
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced_cpp(
	GDALRasterWrapper *p_raster,
	std::vector<size_t> bandIndexes,
	std::string method,
	py::buffer prob,
	std::string filename)
{
	return balanced(
		p_raster,
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
 *
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced_access_cpp(
	GDALRasterWrapper *p_raster,
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
 *
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
balanced_strata_cpp(
	GDALRasterWrapper *p_raster,
	std::vector<size_t> bandIndexes,
	GDALRasterWrapper *p_sraster,
	size_t stratBand,
	std::string method,
	py::buffer prob,
	std::string filename)
{
	return balanced(
		p_raster,
		bandIndexes,
		p_sraster,
		strataBand,
		nullptr,
		"",
		0,
		0,
		method.
		prob,
		filename
	);
}


PYBIND11_MODULE(balanced, m) {
	m.def("balanced_cpp", &balanced_cpp);
	m.def("balanced_access_cpp", &balanced_access_cpp);
	m.def("balanced_strata_cpp", &balanced_strata_cpp);
	m.def("balanced_access_strata_cpp", &balanced)
}
