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

	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using band
		p_mask = maskInfo.second; //pointer to mask
	}

	//allocate spread, balanced, and potentially strata matrices, error check their allocation
	double * p_balanced = (double *) VSIMalloc3(height, width, bandCount * sizeof(double));
	//double * p_balancedCheck = (double *) VSIMalloc3(height, width, bandCount * sizeof(double)); // for testing purposes
	double * p_spread = (double *) VSIMalloc3(height, width, 2 * sizeof(double));
	int * p_strata = (int *)((method == "lcubestratified") ? VSIMalloc3(height, width, sizeof(int)) : nullptr);
	if (!p_balanced || !p_spread || (method == "lcubestratified" && !p_strata)) {
		throw std::runtime_error("unable to allocate raster band to memory");
	}	

	size_t dataPointerIncrement;
	size_t pixelSpace;
	size_t lineSpace;
	if (method == "lpm2_kdtree") {
		dataPointerIncrement = height * width * sizeof(double);
		pixelSpace = 0; //when 0 GDAL sets to appropriate value
		lineSpace = 0; //when 0 GDAL sets to appropriate value
	}
	else {
		dataPointerIncrement = sizeof(double);
		pixelSpace = bandCount * sizeof(double);
		lineSpace = width * bandCount * sizeof(double);
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
			(void *)((size_t)p_balanced + i * dataPointerIncrement),
			width,
			height,
			GDT_Float64,
			pixelSpace,
			lineSpace
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

	//for testing purposes
	/*
	std::memcpy(
		p_balancedCheck,
		p_balanced,
		height * width * bandCount * sizeof(double)
	);
	*/

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

			if (method == "lcubestratified" || method == "lcube") {
				//for lcube methods there are N rows of pixels, and columns equal to the number of bands
				size_t readPixelStart = readIndex * bandCount;
				size_t writePixelStart = writeIndex * bandCount;
				for (size_t j = 0; j < bandCount; j++) {
					double val = p_balanced[readPixelStart + j];
					isNan |= (std::isnan(val) || val == bandNoDataValues[j]);
					if (readIndex != writeIndex) {
						p_balanced[writePixelStart + j] = val;
					}
				}
			}
			else {
				//for lpm2_kdtree method there rows equal to the number of bands, and N columns of pixels
				for (size_t j = 0; j < bandCount; j++) {
					size_t bandStart = j * height * width;
					size_t readPixel = bandStart + readIndex;
					size_t writePixel = bandStart + writeIndex;
					double val = p_balanced[readPixel];
					isNan |= (std::isnan(val) || val == bandNoDataValues[j]);
					if (readIndex != writeIndex) {
						p_balanced[writePixel] = val;
					}
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

	if (method == "lpm2_kdtree") {
		for (size_t band = 1; band < bandCount; band++) {
			void *dst = (void *)((size_t)p_balanced + band * dataPixelCount * sizeof(double));
			void *src = (void *)((size_t)p_balanced + band * height * width * sizeof(double));

			size_t toCopy = dataPixelCount * sizeof(double);
			size_t maxChunkSize = (size_t)src - (size_t)dst;

			while (toCopy > 0) {
				size_t copyAmount = std::min(maxChunkSize, toCopy);
				std::memcpy(
					dst,
					src,
					copyAmount
				);
				dst = (void *)((size_t)dst + copyAmount);
				src = (void *)((size_t)src + copyAmount);
				toCopy -= copyAmount;
			}
		}
	}

	//these for loops are used to ensure the data is copied to buffers correctly
	//they are only to be used for testing.
	/*
	if (method == "lpm2_kdtree") {
		for (size_t band = 0; band < bandCount; band++) {
			size_t start = band * dataPixelCount;
			size_t checkStart = band * height * width;

			size_t checkIndex = 0;
			for (size_t i = 0; i < dataPixelCount; i++) {
				double checkVal = p_balancedCheck[checkStart + checkIndex];
				while (std::isnan(checkVal) || checkVal == bandNoDataValues[band]) {
					checkIndex++;
					checkVal = p_balancedCheck[checkStart + checkIndex];
				}
				double val = p_balanced[start + i];
				if (checkVal != val) {
					std::cout << "[ " << checkIndex << "] = " << checkVal << " NOT EQUAL TO [" << i << "] = " << val << std::endl;
				}
				checkIndex++;
			}
		}
	}
	else {
		for (size_t band = 0; band < bandCount; band++) {

			size_t checkIndex = 0;
			for (size_t i = 0; i < dataPixelCount; i++) {
				double val = p_balanced[i * bandCount + band];
				double checkVal = p_balancedCheck[checkIndex * bandCount + band];
				while(std::isnan(checkVal) || checkVal == bandNoDataValues[band]) {
					checkIndex++;
					checkVal = p_balancedCheck[checkIndex * bandCount + band];
				}
				if (checkVal != val) {
					std::cout << "[ " << checkIndex << "] = " << checkVal << " NOT EQUAL TO [" << i << "] = " << val << std::endl;
				}
				checkIndex++;
			}
		}
	}
	*/

	p_balanced = (double *)VSIRealloc(p_balanced, dataPixelCount * bandCount * sizeof(double));
	p_spread = (double *)VSIRealloc(p_spread, dataPixelCount * 2 * sizeof(double));
	if (!p_balanced || !p_spread) {
		throw std::runtime_error("unable to re allocate rasters");
	}

	double *p_prob;
	std::vector<double> probVect;
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
	Lpm * p_lpm = nullptr;
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
	else if (method == "lpm2_kdtree") {
		p_lpm = new Lpm(
			LpmMethod::lpm2,
			p_prob,
			p_balanced,
			dataPixelCount,
			bandCount,
			eps,
			treeBucketSize,
			treeMethod
		);
		p_lpm->Run();
		indexes = &(p_lpm->sample);
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
	delete p_lpm;

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
