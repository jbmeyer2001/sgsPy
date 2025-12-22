/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#include "utils/raster.h"
#include "utils/helper.h"

#include <condition_variable>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <mkl.h>

namespace sgs {
namespace quantiles {

/* This helper function is used to calculate the quantiles of a
 * raster which is entirely in memory. This is the single
 * precision version of this function, which is used for all 
 * raster data types except double precision floating point 
 * values.
 *
 * Both the quantile probabilities, and the quantiles themselves
 * are double precision, however, so they must be converted
 * before and after execution to and from single precision
 * floating point.
 *
 * Intel's MKL (Math Kernel Library) package is used to calculate the quantiles.
 * The MKL package does not account for nan values, 
 * so the data must be filtered first, leaving only the data 
 * pixel values left.
 *
 * the following MKL VSL (Vector Statistics Library) functions
 * are used to calculate the quantiles:
 * vslsSSCreateTask()
 * vslsSSEditQuantiles()
 * vslsSSCompute()
 * vslSSDeleteTask()
 *
 * @param GDALRasterWrapper *p_raster
 * @param RasterBandMetaData& band
 * @param std::vector<double>& probabilities
 * @param std::vector<double>& quantiles
 */
void calcSPQuantiles(
	raster::GDALRasterWrapper *p_raster,
	helper::RasterBandMetaData& band,
	std::vector<double>& probabilities,
	std::vector<double>& quantiles) 
{
	std::vector<float> spProbabilities(probabilities.size());
	std::vector<float> spQuantiles(quantiles.size());

	//get float vector of probabilities
	for (size_t i = 0; i < probabilities.size(); i++) {
		spProbabilities[i] = static_cast<float>(probabilities[i]);
	}

	//filter out nan values
	size_t fi = 0;
	float nan = static_cast<float>(band.nan);
	size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) * static_cast<size_t>(p_raster->getHeight());
	std::vector<float> filteredData(pixelCount);
	for (size_t i = 0; i < pixelCount; i++) {
		float val = helper::getPixelValueDependingOnType<float>(band.type, band.p_buffer, i);
		bool isNan = std::isnan(val) || val == nan;
		if (!isNan) {
			filteredData[fi] = val;
			fi ++;
		}
	}
	filteredData.resize(fi + 1);

	//define variables to pass to MKL quantiles calculation function
	//https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vslsseditquantiles.html
	VSLSSTaskPtr task;
	int status;
       	MKL_INT quant_order_n = spQuantiles.size();	
	MKL_INT p = 1;
	MKL_INT nparams = filteredData.size();
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;
	float *quant_order = spProbabilities.data();
	float *params = filteredData.data();
	float *quants = spQuantiles.data();

	//calculate quantiles using MKL
	status = vslsSSNewTask(&task, &p, &nparams, &xstorage, params, nullptr, nullptr);
	status = vslsSSEditQuantiles(task, &quant_order_n, quant_order, quants, nullptr, nullptr);
	status = vslsSSCompute(task, VSL_SS_QUANTS, VSL_SS_METHOD_FAST);
	status = vslSSDeleteTask(&task);

	//assign double quantile values
	for(size_t i = 0; i < spQuantiles.size(); i++) {
		quantiles[i] = static_cast<double>(spQuantiles[i]);
	}
}

/* This helper function is used to calculate the quantiles of a
 * raster which is entirely in memory. This is the double
 * precision version of this function.
 *
 * Intel's MKL (Math Kernel Library) package is used to calculate the quantiles.
 * The MKL package does not account for nan values, 
 * so the data must be filtered first, leaving only the data 
 * pixel values left.
 *
 * the following MKL VSL (Vector Statistics Library) functions
 * are used to calculate the quantiles:
 * vsldSSCreateTask()
 * vsldSSEditQuantiles()
 * vsldSSCompute()
 * vslSSDeleteTask()
 *
 * @param GDALRasterWrapper *p_raster
 * @param RasterBandMetaData& band
 * @param std::vector<double>& probabilities
 * @param std::vector<double>& quantiles
 */
void calcDPQuantiles(
	raster::GDALRasterWrapper *p_raster,
	helper::RasterBandMetaData& band,
	std::vector<double>& probabilities,
	std::vector<double>& quantiles) 
{
	//filter out nan values
	size_t fi = 0;
	size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) * static_cast<size_t>(p_raster->getHeight());
	std::vector<double> filteredData(pixelCount);
	for (size_t i = 0; i < pixelCount; i++) {
		double val = helper::getPixelValueDependingOnType<double>(band.type, band.p_buffer, i);
		bool isNan = std::isnan(val) || val == band.nan;
		if (!isNan) {
			filteredData[fi] = val;
			fi ++;
		}
	}
	filteredData.resize(fi + 1);

	//define variables to pass to MKL quantiles calculation function
	//https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vslsseditquantiles.html
	VSLSSTaskPtr task;
	int status;
       	MKL_INT quant_order_n = quantiles.size();	
	MKL_INT p = 1;
	MKL_INT nparams = filteredData.size();
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;
	double *quant_order = probabilities.data();
	double *params = filteredData.data();
	double *quants = quantiles.data();

	//calculate quantiles using MKL
	status = vsldSSNewTask(&task, &p, &nparams, &xstorage, params, nullptr, nullptr);
	status = vsldSSEditQuantiles(task, &quant_order_n, quant_order, quants, nullptr, nullptr);
	status = vsldSSCompute(task, VSL_SS_QUANTS, VSL_SS_METHOD_FAST);
	status = vslSSDeleteTask(&task);
}

/* This helper function is used to calculate the quantiles of a
 * large raster which is more efficient to calculate in batches
 * rather than trying to allocate into memory. This is the single
 * precision version of this function, which is used for all 
 * raster data types except double precision floating point 
 * values.
 *
 * Both the quantile probabilities, and the quantiles themselves
 * are double precision, however, so they must be converted
 * before and after execution to and from single precision
 * floating point.
 *
 * Intel's MKL (Math Kernel Library) package is used to calculate the quantiles.
 * The MKL package does not account for nan values, 
 * so the data must be filtered first, leaving only the data 
 * pixel values left.
 *
 * the following MKL VSL (Vector Statistics Library) functions
 * are used to calculate the quantiles:
 * vslsSSCreateTask()
 * vslsSSEditStreamQuantiles()
 * vslsSSCompute()
 * vslSSDeleteTask()
 *
 * CreateTask is called at the start of execution, and
 * VSLsSSEditStreamQuantiles() is called after the first
 * set of values has been written into the buffer. It was 
 * found that calling this function before resulted in 
 * incorrect values, although the funciton still only
 * has to be called once. The Compute function is called many
 * times, to continuously update the quantiles across batches.
 * The fast calculation method is used, which means afterwards
 * the Compute funciton must be called again with the normal
 * method and an observable count of 0 to get the final quantile
 * values. Quantile streaming algorithms are not exact, since the 
 * entire raster must be in memory to get the most precise
 * possible values. However, the amount of error can be controlled
 * to a specific epsilon (eps) value.
 *
 * The Quantile streaming method is the method introduced by 
 * Zhang et al. and utilized by MKL:
 * https://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf
 * https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-summary-statistics-notes/2021-1/computing-quantiles-with-vsl-ss-method-squants-zw.html
 *
 * Due to the fact that large raster processing in batches is
 * multi-threaded in SGS, a condition variable (cv) is used to 
 * ensure a thread which would be using the calculated quantiles
 * doesn't begin processing or writing any data before the quantiles
 * are finalized. The shared resource is a boolean representing the
 * band which is set to true by this thread once the quantiles are
 * calculated. Any number of processing threads may be waiting on 
 * the quantiles of this band to be calculated, so notify_all is 
 * called on the cv to wake those threads up.
 *
 * @param GDALRasterWrapper *p_raster
 * @param RasterBandMetaData& band
 * @param std::vector<double>& probabilities
 * @param std::vector<double>& quantiles
 * @param std::mutex& mutex
 * @param std::condition_variable& cdv
 * @param bool& calculated
 * @param double eps
 */
void batchCalcSPQuantiles(
	raster::GDALRasterWrapper *p_raster, 
	helper::RasterBandMetaData& band, 
	std::vector<double>& probabilities,
	std::vector<double>& quantiles,
	std::mutex& mutex,
	std::condition_variable& cv,
	bool& calculated,
	double eps) 
{
	std::vector<float> spProbabilities(probabilities.size());
	std::vector<float> spQuantiles(quantiles.size());

	//get float vector of probabilities
	for (size_t i = 0; i < probabilities.size(); i++) {
		spProbabilities[i] = static_cast<float>(probabilities[i]);
	}

	int xBlocks = (p_raster->getWidth() + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (p_raster->getHeight() + band.yBlockSize - 1) / band.yBlockSize;
	
	float nan = static_cast<float>(band.nan);
	float spEps = static_cast<float>(eps);
	void *p_buffer = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);
	float *p_filtered = reinterpret_cast<float *>(VSIMalloc3(band.xBlockSize, band.yBlockSize, sizeof(float)));

	//define variables to pass to MKL quantiles calculation function
	//https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vslsseditstreamquantiles.html
	VSLSSTaskPtr task;
	int status;
       	MKL_INT quant_order_n = spProbabilities.size();	
	MKL_INT p = 1;
	MKL_INT n = band.xBlockSize * band.yBlockSize;
	MKL_INT nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;

	float *quant_order = reinterpret_cast<float *>(spProbabilities.data());
	float *quants = reinterpret_cast<float *>(spQuantiles.data());

	status = vslsSSNewTask(&task, &p, &n, &xstorage, p_filtered, nullptr, nullptr);

	//read and compute raster by blocks
	bool calledEditStreamQuantiles = false;
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_mutex->lock();
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			CPLErr err = band.p_band->ReadBlock(xBlock, yBlock, p_buffer);
			band.p_mutex->unlock();
			if (err) {
				throw std::runtime_error("error reading block from band.");
			}

			int fi = 0;
			for (int y = 0; y < yValid; y++) {
				int index = y * band.xBlockSize;
				for (int x = 0; x < xValid; x++) {
					float val = helper::getPixelValueDependingOnType<float>(band.type, p_buffer, index);
					bool isNan = std::isnan(val) || val == nan;
					if (!isNan) {
						p_filtered[fi] = val;
						fi++;
					}
					index++;
				}
			}

			if (fi == 0) {
				continue;
			}

			n = fi;
			if (!calledEditStreamQuantiles) {
				//for some reason this has issues if called before the first band of data is loaded,
				//although it only has to be called once.
				status = vslsSSEditStreamQuantiles(task, &quant_order_n, quant_order, quants, &nparams, &spEps);
				calledEditStreamQuantiles = true;
			}
			status = vslsSSCompute(task, VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZW_FAST);
		}
	}

	//use VSL_SS_METHOD_SQUANTS_ZW (not VSL_SS_METHOD_SQUANTS_ZW_FAST) to get final estimates
	n = 0;
	vslsSSCompute(task, VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZW);	
	status = vslSSDeleteTask(&task);

	//assign double quantile values
	for(size_t i = 0; i < spQuantiles.size(); i++) {
		quantiles[i] = static_cast<double>(spQuantiles[i]);
	}

	VSIFree(p_buffer);
	VSIFree(p_filtered);

	mutex.lock();
	calculated = true;
	mutex.unlock();
	cv.notify_all();
}


/* This helper function is used to calculate the quantiles of a
 * large raster which is more efficient to calculate in batches
 * rather than trying to allocate into memory. This is the double
 * precision version of this function.
 *
 * Intel's MKL (Math Kernel Library) package is used to calculate the quantiles.
 * The MKL package does not account for nan values, 
 * so the data must be filtered first, leaving only the data 
 * pixel values left.
 *
 * the following MKL VSL (Vector Statistics Library) functions
 * are used to calculate the quantiles:
 * vsldSSCreateTask()
 * vsldSSEditStreamQuantiles()
 * vsldSSCompute()
 * vslSSDeleteTask()
 *
 * CreateTask is called at the start of execution, and
 * VSLsSSEditStreamQuantiles() is called after the first
 * set of values has been written into the buffer. It was 
 * found that calling this function before resulted in 
 * incorrect values, although the funciton still only
 * has to be called once. The Compute function is called many
 * times, to continuously update the quantiles across batches.
 * The fast calculation method is used, which means afterwards
 * the Compute funciton must be called again with the normal
 * method and an observable count of 0 to get the final quantile
 * values. Quantile streaming algorithms are not exact, since the 
 * entire raster must be in memory to get the most precise
 * possible values. However, the amount of error can be controlled
 * to a specific epsilon (eps) value.
 *
 * The Quantile streaming method is the method introduced by 
 * Zhang et al. and utilized by MKL:
 * https://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf
 * https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-summary-statistics-notes/2021-1/computing-quantiles-with-vsl-ss-method-squants-zw.html
 *
 * Due to the fact that large raster processing in batches is
 * multi-threaded in SGS, a condition variable (cv) is used to 
 * ensure a thread which would be using the calculated quantiles
 * doesn't begin processing or writing any data before the quantiles
 * are finalized. The shared resource is a boolean representing the
 * band which is set to true by this thread once the quantiles are
 * calculated. Any number of processing threads may be waiting on 
 * the quantiles of this band to be calculated, so notify_all is 
 * called on the cv to wake those threads up.
 *
 * @param GDALRasterWrapper *p_raster
 * @param RasterBandMetaData& band
 * @param std::vector<double>& probabilities
 * @param std::vector<double>& quantiles
 * @param std::mutex& mutex
 * @param std::condition_variable& cdv
 * @param bool& calculated
 * @param double eps
 */
void batchCalcDPQuantiles(
	raster::GDALRasterWrapper *p_raster, 
	helper::RasterBandMetaData& band, 
	std::vector<double>& probabilities,
	std::vector<double>& quantiles,
	std::mutex& mutex,
	std::condition_variable& cv,
	bool& calculated,
	double eps) 
{
	int xBlocks = (p_raster->getWidth() + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (p_raster->getHeight() + band.yBlockSize - 1) / band.yBlockSize;
	
	void *p_buffer = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);
	double *p_filtered = reinterpret_cast<double *>(VSIMalloc3(band.xBlockSize, band.yBlockSize, sizeof(double)));

	//define variables to pass to MKL quantiles calculation function
	//https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vslsseditstreamquantiles.html
	VSLSSTaskPtr task;
	int status;
       	MKL_INT quant_order_n = probabilities.size();	
	MKL_INT p = 1;
	MKL_INT n = band.xBlockSize * band.yBlockSize;
	MKL_INT nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;

	double *quant_order = reinterpret_cast<double *>(probabilities.data());
	double *quants = reinterpret_cast<double *>(quantiles.data());

	vsldSSNewTask(&task, &p, &n, &xstorage, p_filtered, 0, 0);

	//read and compute raster by blocks
	bool calledEditStreamQuantiles = false;
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_mutex->lock();
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			CPLErr err = band.p_band->ReadBlock(xBlock, yBlock, p_buffer);
			band.p_mutex->unlock();
			if (err) {
				throw std::runtime_error("error reading block from band.");
			}

			int fi = 0;
			for (int y = 0; y < yValid; y++) {
				int index = y * band.xBlockSize;
				for (int x = 0; x < xValid; x++) {
					double val = helper::getPixelValueDependingOnType<double>(band.type, p_buffer, index);
					bool isNan = std::isnan(val) || val == band.nan;
					if (!isNan) {
						p_filtered[fi] = val;
						fi++;
					}
					index++;
				}
			}

			if (fi == 0) {
				continue;
			}

			n = fi;
			if (!calledEditStreamQuantiles) {
				//for some reason this has issues if called before the first band of data is loaded,
				//although it only has to be called once.
				status = vsldSSEditStreamQuantiles(task, &quant_order_n, quant_order, quants, &nparams, &eps);
				calledEditStreamQuantiles = true;
			}
			status = vsldSSCompute(task, VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZW_FAST);
		}
	}

	//use VSL_SS_METHOD_SQUANTS_ZW (not VSL_SS_METHOD_SQUANTS_ZW_FAST) to get final estimates
	n = 0;
	vsldSSCompute(task, VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZW);	
	status = vslSSDeleteTask(&task);

	mutex.lock();
	calculated = true;
	mutex.unlock();
	cv.notify_all();

	VSIFree(p_buffer);
	VSIFree(p_filtered);
}

/**
 * This is a helper function for processing a pixel of data
 * when a mapped stratification is being created.
 *
 * First, the value is read in as a double, and it is determined
 * whether the pixel is a nan pixel or not. The mapNan boolean 
 * is updated in addition to the isNan boolean, to ensure that
 * if one band within the raster is nan at a certain pixel
 * then the mapped raster (but not necessarily all output rasters)
 * is also nan at that pixel.
 *
 * Then, if it isn't a nan pixel the lower bound of the value 
 * within the vector of break values is found. For example,
 * if the value was 3 and the breaks vector was [2, 4, 6], the
 * lower bound would be 1, which is the index of 4, the first
 * value larger than 3 in the breaks vector. This lower bound
 * is the strata. This strata (or the nan value) is then written 
 * with the appropriate type to the strat raster band.
 *
 * @param size_t index
 * @param RasterBandMetaData& dataBand
 * @param void *p_dataBuffer
 * @param RasterBandMetaData& stratBand
 * @param void *p_stratBuffer
 * @param std::vector<double>& quantiles
 * @param size_t multiplier
 * @param bool& mapNan
 * @param size_t &mapStrat
 */
inline void processMapPixel(
	size_t index,
	helper::RasterBandMetaData& dataBand,
	void * p_dataBuffer,
	helper::RasterBandMetaData& stratBand,
	void * p_stratBuffer,
	std::vector<double>& quantiles,
	size_t multiplier,
	bool& mapNan,
	size_t& mapStrat)
{
	double val = helper::getPixelValueDependingOnType<double>(dataBand.type, p_dataBuffer, index);
	bool isNan = std::isnan(val) || (double)val == dataBand.nan;
	mapNan |= isNan;
		
	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(quantiles.begin(), quantiles.end(), val);
		strat = (it == quantiles.end()) ? quantiles.size() : std::distance(quantiles.begin(), it);
	}

	helper::setStrataPixelDependingOnType(stratBand.type, p_stratBuffer, index, isNan, strat);

	if (!mapNan) {
		mapStrat += strat * multiplier;
	}
}	

/**
 * This is a helper function for processing a pixel of data.
 *
 * First, the value is read in as a double, and it is determined
 * whether the pixel is a nan pixel or not.
 *
 * Then, if it isn't a nan pixel the lower bound of the value 
 * within the vector of break values is found. For example,
 * if the value was 3 and the breaks vector was [2, 4, 6], the
 * lower bound would be 1, which is the index of 4, the first
 * value larger than 3 in the breaks vector. This lower bound
 * is the strata. This strata (or the nan value) is then written 
 * with the appropriate type to the strat raster band.
 *
 * @param size_t index
 * @param void *p_data
 * @param RasterBandMetaData *p_dataBand
 * @param void *p_strat
 * @param RasterBandMetaData *p_stratBand
 * @param std::vector<double>& quantiles
 */
inline void
processPixel(
	size_t index,
	void *p_data,
	helper::RasterBandMetaData *p_dataBand,
	void *p_strat,
	helper::RasterBandMetaData *p_stratBand,
	std::vector<double>& quantiles)
{
	double val = helper::getPixelValueDependingOnType<double>(p_dataBand->type, p_data, index);
	bool isNan = std::isnan(val) || val == p_dataBand->nan;

	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(quantiles.begin(), quantiles.end(), val);
		strat = (it == quantiles.end()) ? quantiles.size() : std::distance(quantiles.begin(), it);
	}

	helper::setStrataPixelDependingOnType(p_stratBand->type, p_strat, index, isNan, strat);
}

/**
 * This function stratifies a given raster using user-defined probabilities.
 * The probabilities (quantiles) are provided as a vector of doubles mapped
 * to a band index.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map variable to combine the stratification of
 * the given raster bands.
 *
 * The function can be thought of in four different sections: the setup,
 * the quantiles calculation, the processing, and the finish/return. During
 * the setup, metadata is acquired for the input raster, and an output dataset
 * is created which depends on user-given parameters and the input raster.
 * During the quantiles calculation, Intel's Math Kernel Library (MKL) is
 * used to calculate quantiles either with the entire raster in memory
 * or in batches. During the processing step the raster is iterated through,
 * either by blocks or with the entire raster in memory, the strata are determined
 * for each pixel and then written to the output dataset. During the finish/return 
 * step, a GDALRasterWrapper object is created using the output dataset.
 *
 *
 * SETUP:
 * the data structures holding metadata are initialized and it is determined
 * whether the raster is a virtual raster or not, and if it is a virtual raster
 * whether it is fully in-memory or whether it must be stored on disk.
 *
 * If the user provides an output filename, the dataset will not be a virtual
 * dataset instead it will be associated with the filename. If the user does
 * not provide an output filename then a virtual dataset driver will be used.
 * In the case of a large raster (whether or not the raster is large enough for
 * this is calculated and passed by hte Python side of the application), the dataset
 * will be VRT. If the package is comfortable fitting the entire raster in memory,
 * an in-memory dataset will be used.
 *
 * The input raster bands are iterated through, metadata is stored on them,
 * and bands are created for the output dataset. In the case of a VRT dataset,
 * each band is a complete dataset itself which must be added after it has
 * been written to. In the case of a MEM dataset, the bands must be acquired
 * from the input raster. Both MEM and VRT support the addBand function,
 * and support bands of different types. Thus, the bands are added dynamically
 * while iterating through input raster bands. Non virtual formats require the
 * data types to be known at dataset initialization and don't support the
 * addBand function, so the dataset must be created after iterating through
 * the input bands.
 *
 *
 * QUANTILES CALCULATION:
 * There are four different functions which do the quantiles calculation.
 * Different functions are used depending on whether the raster should be
 * batch processed, or whether it is entirely in memory. Further, there
 * are different functions for single precision (float) or double precision
 * (double) floating point values. Intels Math Kernel Library (MKL) is used
 * to calculate the quantiles.
 *
 * For single precision vs. double precision, the only case where double
 * precision would be used is if the raster data type is double. the
 * functions within the MKL library are slightly different for double 
 * and single precision. Further, the input and output vectors representing
 * quantile probabilities and the quantiles themselves are vectors of
 * double, so in the case where single precision is used, intermediate
 * versions of the data must be written to and read from at the beginning
 * and end of the function. The reasing why double precision is not used
 * in every case, is because the same amount of values of single precision
 * take up half the ammount of memory as their double precision counterparts.
 * As a result, the memory usage due to intermediate stored values should
 * be significantly less for single precision floating point values.
 *
 * For batch processing vs non-batch processing, slightly different algorithms
 * are used. The batch processing algorithm is a quantile streaming algorithm
 * within MKL which creates a fast and accurate approximation of the quantiles
 * without requiring the whole raster to be in memory at once. More information
 * is contained in the documentation for those particular functions.
 *
 * Information on the quantile streaming algorithm can be found here:
 *  - https://www.cs.unc.edu/~weiwang/paper/SSDBM07_2.pdf
 *  - https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-summary-statistics-notes/2021-1/computing-quantiles-with-vsl-ss-method-squants-zw.html
 *
 * PROCESSING:
 * the processing section iterated through ever pixel in every input band,
 * and calculates/writes the strata to the corresponding output band.
 *
 * There are four different cases dealing with whether or not the entire
 * raster band is allocated in memory (the largeRaster bariable is false),
 * and whether or not the values of each band should be mapped to an extra
 * output raster band.
 *
 * If the raster is large, it is porcessed in blocks and split into
 * groups of blocks to be processed by multiple threads. If the raster 
 * bands are in-memory, the entire raster is processed at once by a 
 * single thread. The mapped rasters store information on an extra
 * output raster band, the output values of which are determined
 * as a function of all other output raster bands. The multipliers
 * vector stores the information for this.
 *
 * For the large rasters, the processing starts out by splitting the
 * raster into chunks depending on the number of threads. A thread is
 * then ctreated for each chunk. Within each thread, the blocks within 
 * it's designated chunk are iterated through and first read from the
 * input bands, processed, then written to the output bands. In the
 * case of a mapped raster all of the bands are iterated through
 * alongside achother so that the intermediate mapping calculations
 * don't have to be written then read again. In the case of a non
 * mapped raster, each band is processed sequentially. 
 *
 * CLEANUP:
 * If the output dataset is a VRT dataset, the dataset which represent
 * its bands (which have not yet been added as bands) must be added as
 * bands now that they are populated with data and are thus allowed
 * to be added.
 *
 * If the dataset output bands are fully in memory, they are moved to a
 * vector from their metadata objects to be passed as a parameter to the
 * GDALRasterWrapper constructor (or not if the bands aren't in memory). 
 * This GDALRasterWrapper is then returned.
 * 
 * @param GDALRasterWrapper *p_raster
 * @param std::map<int, std::vector<double>> userProbabilities
 * @param bool map 
 * @param std::string filename
 * @param std::string tempFolder,
 * @param bool largeRaster,
 * @param int threadCount
 * @param std::map<std::string, std::string> driverOptions,
 * @param double eps
 * @returns GDALRasterWrapper *pointer to newly created stratified raster
 */
raster::GDALRasterWrapper *quantiles(
	raster::GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userProbabilites,
	bool map,
	std::string filename,
	std::string tempFolder,
	bool largeRaster,
	int threadCount,
	std::map<std::string, std::string> driverOptions,
	double eps) 
{
	GDALAllRegister();

	int height = p_raster->getHeight();
	int width = p_raster->getWidth();
	double *geotransform = p_raster->getGeotransform();
	std::string projection = std::string(p_raster->getDataset()->GetProjectionRef());

	GDALDataset *p_dataset = nullptr;

	int bandCount = userProbabilites.size();
	std::vector<std::vector<double>> probabilities;
	std::vector<std::string> bandNames = p_raster->getBands();

	std::vector<helper::RasterBandMetaData> dataBands(bandCount);
	std::vector<helper::RasterBandMetaData> stratBands(bandCount + map);
	std::vector<helper::VRTBandDatasetInfo> VRTBandInfo;

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";

	std::mutex dataBandMutex;
	std::mutex stratBandMutex;
	std::vector<std::mutex> stratBandMutexes(isVRTDataset * (bandCount + map));

	std::string driver;
	if (isMEMDataset || isVRTDataset) {
		driver = isMEMDataset ? "MEM" : "VRT";
	       	p_dataset = helper::createVirtualDataset(driver, width, height, geotransform, projection);	
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();

		if (extension == ".tif") {
			driver = "Gtiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now");
		}
	}

	//allocate, read, and initialize raster data and breaks information
	GDALDataType stratPixelType = GDT_Int8;
	size_t stratPixelSize = 1;
	size_t band = 0;
	for (auto const& [key, val] : userProbabilites) {
		helper::RasterBandMetaData *p_dataBand = &dataBands[band];
		helper::RasterBandMetaData *p_stratBand = &stratBands[band];

		//get and store metadata from input raster band
		GDALRasterBand *p_band = p_raster->getRasterBand(key);
		p_dataBand->p_band = p_band;
		p_dataBand->type = p_raster->getRasterBandType(key);
		p_dataBand->size = p_raster->getRasterBandTypeSize(key);
		p_dataBand->p_buffer = largeRaster ? nullptr : p_raster->getRasterBandBuffer(key);
		p_dataBand->nan = p_band->GetNoDataValue();
		p_dataBand->p_mutex = &dataBandMutex;
		p_band->GetBlockSize(&p_dataBand->xBlockSize, &p_dataBand->yBlockSize);

		//add probabilities
		probabilities.push_back(val);	

		//update metadata of new strat raster
		size_t maxStrata = val.size() + 1;
		helper::setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_" + bandNames[key];
		p_stratBand->xBlockSize = map ? dataBands[0].xBlockSize : p_dataBand->xBlockSize;
		p_stratBand->yBlockSize = map ? dataBands[0].yBlockSize : p_dataBand->yBlockSize;
		p_stratBand->p_mutex = isVRTDataset ? &stratBandMutexes[band] : &stratBandMutex;

		//update dataset with new band information
		if (isMEMDataset) {
			helper::addBandToMEMDataset(p_dataset, *p_stratBand);
		}
		else if (isVRTDataset) {
			helper::createVRTBandDataset(p_dataset, *p_stratBand, tempFolder, std::to_string(key), VRTBandInfo, driverOptions);
		}
		else {
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}

		band++;
	}

	//set multipliers if mapped stratification
	std::vector<size_t> multipliers(probabilities.size(), 1);
	if (map) {
		//determine the stratification band index multiplier for the mapped band
		for (int i = 0; i < bandCount - 1; i++) {
			multipliers[i + 1] = multipliers[i] * (probabilities[i].size() + 1);
		}

		//update info of new strat raster band map
		helper::RasterBandMetaData *p_stratBand = &stratBands.back();
		size_t maxStrata = multipliers.back() * (probabilities.back().size() + 1);
		helper::setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_map";
		p_stratBand->xBlockSize = dataBands[0].xBlockSize;
		p_stratBand->yBlockSize = dataBands[0].yBlockSize;
		p_stratBand->p_mutex = isVRTDataset ? &stratBandMutexes.back() : &stratBandMutex;

		//update dataset with band information
		if (isMEMDataset) {
			helper::addBandToMEMDataset(p_dataset, *p_stratBand);
		}
		else if (isVRTDataset) {
			helper::createVRTBandDataset(p_dataset, *p_stratBand, tempFolder, "map", VRTBandInfo, driverOptions);
		}
		else {
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}
	}

	//create full non-virtual dataset now that we have all required band information
	if (!isMEMDataset && !isVRTDataset) {
		bool useTiles = stratBands[0].xBlockSize != width &&
				stratBands[0].yBlockSize != height;

		for (size_t band = 0; band < stratBands.size(); band++) {
			stratBands[band].size = stratPixelSize;
			stratBands[band].type = stratPixelType;
			stratBands[band].p_buffer = !largeRaster ? VSIMalloc3(height, width, stratPixelSize) : nullptr;
		}

		p_dataset = helper::createDataset(
			filename,
			driver,
			width,
			height,
			geotransform,
			projection,
			stratBands.data(),
			stratBands.size(),
			useTiles,
			driverOptions
		);
	}

	std::vector<std::vector<double>> quantiles(probabilities.size());
	if (largeRaster) {
		pybind11::gil_scoped_acquire acquire;
		boost::asio::thread_pool pool(threadCount); 
	
		//initialize synchronization variables
		std::vector<std::mutex> mutexes(map ? 1 : bandCount);
		std::vector<std::condition_variable> cvs(map ? 1 : bandCount);
		bool *quantilesCalculated = reinterpret_cast<bool *>(VSIMalloc2(bandCount, sizeof(bool)));

		//call batch processing quantiles function depending on data type
		for (int i = 0; i < bandCount; i++) {
			helper::RasterBandMetaData band = dataBands[i];
			quantiles[i].resize(probabilities[i].size());
			quantilesCalculated[i] = false;
			if (band.type != GDT_Float64) {
				boost::asio::post(pool, std::bind(batchCalcSPQuantiles,
					p_raster,
					std::ref(band),
					std::ref(probabilities[i]),
					std::ref(quantiles[i]),
					std::ref(mutexes[map ? 0 : i]),
					std::ref(cvs[map ? 0 : i]),
					std::ref(quantilesCalculated[i]),
					static_cast<double>(eps)
				));
			}
			else {
				boost::asio::post(pool, std::bind(batchCalcDPQuantiles,
					p_raster,
					std::ref(band),
					std::ref(probabilities[i]),
					std::ref(quantiles[i]),
					std::ref(mutexes[i]),
					std::ref(cvs[i]),
					std::ref(quantilesCalculated[i]),
					static_cast<double>(eps)
				));
			}
		}

		//iterate through all pixels and update the stratified raster bands
		if (map) {

			//NOTE: all processing threads in the case of a mapped raster must wait
			//until quantiles have been calculated for every band.
			std::unique_lock lock(mutexes[0]);
			bool allBandsQuantilesCalculated = true;
			for (int i = 0; i < bandCount; i++) {
				allBandsQuantilesCalculated &= quantilesCalculated[i];
			}

			while (!allBandsQuantilesCalculated) {
				cvs[0].wait(lock);
				allBandsQuantilesCalculated = true;
				for (int i = 0; i < bandCount; i++) {
					allBandsQuantilesCalculated &= quantilesCalculated[i];
				}
			}

			//use the first raster band to determine block size
			int xBlockSize = dataBands[0].xBlockSize;
			int yBlockSize = dataBands[0].yBlockSize;

			int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
			int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;
			int chunkSize = yBlocks / threadCount;

			for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
				int yBlockEnd = std::min(yBlockStart + chunkSize, yBlocks);
			
				boost::asio::post(pool, [
					bandCount,
					xBlockSize, 
					yBlockSize, 
					yBlockStart, 
					yBlockEnd, 
					xBlocks, 
					&dataBands, 
					&stratBands, 
					&quantiles,
					&multipliers
				] {			
					std::vector<void *> dataBuffers(dataBands.size());
					std::vector<void *> stratBuffers(stratBands.size());
					for (size_t band = 0; band < dataBuffers.size(); band++) {
						dataBuffers[band] = VSIMalloc3(xBlockSize, yBlockSize, dataBands[band].size);
						stratBuffers[band] = VSIMalloc3(xBlockSize, yBlockSize, stratBands[band].size);
					}
					stratBuffers.back() = VSIMalloc3(xBlockSize, yBlockSize, stratBands.back().size);

					for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
						for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
							int xValid, yValid;
							dataBands[0].p_mutex->lock();
							dataBands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
							dataBands[0].p_mutex->unlock();

							//read raster band data into band buffers
							for (size_t band = 0; band < static_cast<size_t>(bandCount); band++) {
								helper::rasterBandIO(
									dataBands[band], 
									dataBuffers[band], 
									xBlockSize, 
									yBlockSize, 
									xBlock, 
									yBlock, 
									xValid, 
									yValid, 
									true //read = true
								);		
							}

							//process blocked band data
							for (int y = 0; y < yValid; y++) {
								size_t index = static_cast<size_t>(y * xBlockSize);
								for (int x = 0; x < xValid; x++) {
									bool mapNan = false;
									size_t mapStrat = 0;
									
									for (size_t band = 0; band < static_cast<size_t>(bandCount); band++) {
										processMapPixel(
											index, 
											dataBands[band], 
											dataBuffers[band], 
											stratBands[band], 
											stratBuffers[band], 
											quantiles[band],
											multipliers[band],
											mapNan,
											mapStrat
										);
									}
								
									helper::setStrataPixelDependingOnType(
										stratBands.back().type,
										stratBuffers.back(),
										index,
										mapNan,
										mapStrat
									);

									index++;
								}
							}
					
							//write strat band data
							for (size_t band = 0; band <= static_cast<size_t>(bandCount); band++) {
								helper::rasterBandIO(
									stratBands[band],
									stratBuffers[band],
									xBlockSize,
									yBlockSize,
									xBlock,
									yBlock,
									xValid,
									yValid,
									false //read = false
								);
							}
						}
					}

					for (size_t band = 0; band < dataBuffers.size(); band++) {
						VSIFree(dataBuffers[band]);
						VSIFree(stratBuffers[band]);
					}
					VSIFree(stratBuffers.back());
				});
			}	
		}
		else {	
			for (size_t band = 0; band < static_cast<size_t>(bandCount); band++) {
				helper::RasterBandMetaData* p_dataBand = &dataBands[band];
				helper::RasterBandMetaData* p_stratBand = &stratBands[band];

				int xBlockSize = p_dataBand->xBlockSize;
				int yBlockSize = p_dataBand->yBlockSize;
					
				int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
				int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;			
				int chunkSize = yBlocks / threadCount;
				
				for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
					int yBlockEnd = std::min(yBlocks, yBlockStart + chunkSize);
					std::vector<double> *p_quantiles = &quantiles[band];
					std::mutex *p_mutex = &mutexes[band];
					std::condition_variable *p_cv = &cvs[band];
					bool *p_calculated = &quantilesCalculated[band];
					
					boost::asio::post(pool, [
						xBlockSize, 
						yBlockSize, 
						yBlockStart, 
						yBlockEnd, 
						xBlocks, 
						p_dataBand, 
						p_stratBand, 
						p_quantiles,
						p_mutex,
						p_cv,
						p_calculated
					] {
						std::unique_lock lock(*p_mutex);

						if (!(*p_calculated)) {
							p_cv->wait(lock);
						}

						void *p_data = VSIMalloc3(xBlockSize, yBlockSize, p_dataBand->size);
						void *p_strat = VSIMalloc3(xBlockSize, yBlockSize, p_stratBand->size);

						int xValid, yValid;
						for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
							for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
								p_dataBand->p_mutex->lock();
								p_dataBand->p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
								p_dataBand->p_mutex->unlock();
								
								//read block into memory
								rasterBandIO(
									*p_dataBand,
									p_data,
									xBlockSize,
									yBlockSize,
									xBlock,
									yBlock,
									xValid,
									yValid,
									true //read = true
								);

								//process block
								for (int y = 0; y < yValid; y++) {
									size_t index = static_cast<size_t>(y * xBlockSize);
									for (int x = 0; x < xValid; x++) {
										processPixel(index, p_data, p_dataBand, p_strat, p_stratBand, *p_quantiles);
										index++;
									}
								}
								
								//write resulting stratifications to disk
								helper::rasterBandIO(
									*p_stratBand,
									p_strat,
									xBlockSize,
									yBlockSize,
									xBlock,
									yBlock,
									xValid,
									yValid,
									false //read = false
								);
							}
						}
						VSIFree(p_data);
						VSIFree(p_strat);
					});
				}
			}
		}

		pool.join();
		VSIFree(quantilesCalculated);
		pybind11::gil_scoped_release release;
	}
	else {
		//call quantiles calculation fuction depending on type
		for (int i = 0; i < bandCount; i++) {
			helper::RasterBandMetaData band = dataBands[i];
			quantiles[i].resize(probabilities[i].size());
			(band.type != GDT_Float64) ?
				calcSPQuantiles(p_raster, band, probabilities[i], quantiles[i]) :
				calcDPQuantiles(p_raster, band, probabilities[i], quantiles[i]);
		}

		size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) * static_cast<size_t>(p_raster->getHeight());
		if (map) {
			for (size_t index = 0; index < pixelCount; index++) {
				bool mapNan = false;
				size_t mapStrat = 0;
									
				for (size_t band = 0; band < static_cast<size_t>(bandCount); band++) {
					processMapPixel(
						index, 
						dataBands[band], 
						dataBands[band].p_buffer, 
						stratBands[band], 
						stratBands[band].p_buffer, 
						quantiles[band],
						multipliers[band],
						mapNan,
						mapStrat
					);
				}

				helper::setStrataPixelDependingOnType(
					stratBands.back().type,
					stratBands.back().p_buffer,
					index,
					mapNan,
					mapStrat
				);
			}
		}
		else {
			for (size_t band = 0; band < static_cast<size_t>(bandCount); band++) {
				for (size_t i = 0; i < pixelCount; i++) {
					processPixel(
						i,
					       	dataBands[band].p_buffer,	
						&dataBands[band], 
						stratBands[band].p_buffer,
						&stratBands[band],
						quantiles[band]
					);
				}
			}
		}

		//if non-virtual dataset then write in-memory output bands to disk
		if (!isVRTDataset && !isMEMDataset) {
			CPLErr err;
			for (size_t band = 0; band < stratBands.size(); band++) {
				err = stratBands[band].p_band->RasterIO(
					GF_Write,
					0,
					0,
					width,
					height,
					stratBands[band].p_buffer,
					width,
					height, 
					stratBands[band].type,
					0,
					0
				);
				if (err) {
					throw std::runtime_error("error writing band to file.");
				}
			}
		}

	}
		
	//close and add all of the VRT sub datasets as bands
	if (isVRTDataset) {
		for (size_t band = 0; band < VRTBandInfo.size(); band++) {
			GDALClose(VRTBandInfo[band].p_dataset);
			helper::addBandToVRTDataset(p_dataset, stratBands[band], VRTBandInfo[band]);
		}
	}

	//if bands are in memory, populate a vector of just bands to use in GDALRasterWrapper creation
	std::vector<void *> buffers(stratBands.size());
	if (!largeRaster) {
		for (size_t band = 0; band < stratBands.size(); band++) {
			buffers[band] = stratBands[band].p_buffer;
		}
	}

	return largeRaster ? 
		new raster::GDALRasterWrapper(p_dataset) :
		new raster::GDALRasterWrapper(p_dataset, buffers);
}

} //namespace quantiles
} //namespace sgs
