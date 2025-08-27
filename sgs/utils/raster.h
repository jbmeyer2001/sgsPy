/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for raster operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#pragma once

#include <filesystem>
#include <iostream>

#include <gdal_priv.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//used as cutoff for max band allowed in memory
#define GIGABYTE 1073741824

namespace py = pybind11;
using namespace pybind11::literals;

/**
 * Wrapper class for a GDAL dataset containing a raster image.
 *
 * This class provides getter methods for important raster data and metadata, 
 * as well as a way to access the raster as an array. 
 *
 * Memory is controlled using either a smart pointer (in the case of the dataset)
 * or CPLMalloc and CPLFree (in the case of the raster image). The function
 * GDALDataset::RasterIO() function should do the windowing for us behind the
 * scenes, meaning we don't have to worry about physical memory size when
 * allocating for large raster images.
 *
 * The raster buffer contains all of the raster bands and can be indexed
 * in Python with: [band][y][x] 
 * and in C++ with: [band * size * width * height + y * size * width + x * size] 
 *
 * The buffer is exposed to the Python side of the application using
 * a py::buffer, and it's exposed to the C++ side of the application using
 * a void *. The expectation is that this void pointer will be cast 
 * to another data type pointer as required.
 */
class GDALRasterWrapper {
	private:
	GDALDatasetUniquePtr p_dataset;

	std::vector<void *> rasterBandPointers;
	std::vector<bool> rasterBandRead;

	std::vector<void *> displayRasterBandPointers;
	std::vector<bool> displayRasterBandRead;
	int displayRasterWidth = -1;
	int displayRasterHeight = -1;

	double geotransform[6];
	/**
	 * Internal function used to read raster band data.
	 * 
	 * the band int is used to select the band from the dataset.
	 *
	 * width and height are used to allocate memory and are passed into
	 * the GDALRasterBand::RasterIO function. The band
	 * may be downsampled depending on the values of width and height.
	 *
	 * If this function completes successfully, the void * in either
	 * rasterBandPointers, or displayRasterBandPointers (in the case
	 * where height/width are not the same as the full image) which
	 * corresponds to the band index given will be the allocated
	 * data buffer to the raster band.
	 *
	 * @param int width
	 * @param int height
	 * @param int band zero-indexed
	 * @throws std::runtime_error if unable to read raster band
	 */
	void readRasterBand(int width, int height, int band);

	/**
	 * Internal function which returns a pybuffer of the raster band, using
	 * the type specified. The band must have already been allocated
	 * using allocateRasterBand(). 
	 *
	 * @param size_t size of data type for one pixel
	 * @param void *p_raster pointer to allocated raster
	 * @param int width 
	 * @param int height
	 * @returns py::buffer memoryview object of the data
	 */
	template <typename T> 
	py::buffer getBuffer(size_t size, void *p_raster, int width, int height);
	
	public:
	/**
	 * Constructor for GDALRasterWrapper class.
	 * Creates GDALDataset using drivers and given file, then calls
	 * createFromDataset() passing the created object.
	 *
	 * @param filename as std::string
	 * @throws std::runtime_error if dataset is not initialized
	 * @throws std::runtime_error if unable to get geotransform
	 */
	GDALRasterWrapper(std::string filename);

	/**
	 * Constructor for GDALRasterWrapper class if an in-memory dataset
	 * has already been created. The bands (already read and allocated)
	 * are passed as the second parameter.
	 * Calls createFromDataset() passing p_dataset parameter, and
	 * sets internal raster band parameters accordingly.
	 *
	 * @param GDALDataset *p_dataset GDAL raster dataset
	 * @param std::vector<void *> raster bands
	 * @throws std::runtime_error if unable to get geotransform
	 */
	GDALRasterWrapper(GDALDataset *p_dataset, std::vector<void *> bands);

	/**
	 * Constructor for GDALRasterWrapper class using just a GDAL 
	 * dataset pointer, by calling createFromDataset().
	 *
	 * @param GDALDataset *p_dataset GDAL raster dataset
	 */
	GDALRasterWrapper(GDALDataset *p_dataset);

	/**
	 * Constructor for GDALRasterWrapper class. This method creates
	 * a GDAL dataset either in memory using the raster bands specified.
	 *
	 * @param std::vector<void *> rasterBands the bands which will be part of the raster
	 * @param GDALDataType type the pixel data type
	 * @param double *geotransform of new raster
	 */
	GDALRasterWrapper(
		std::vector<void *> rasterBands, 
		std::vector<std::string> rasterBandNames,
		std::vector<GDALDataType> rasterBandTypes,
		int width, 
		int height, 
		double *geotransform,
		std::string projection
	);

	/**
	 * Populates/constructs the GDALRasterWrapper object using a raster
	 * dataset pointer. Called by some of the GDALRasterWrapper 
	 * constructors as a helper function.
	 *
	 * @param GDALDataset *GDAL raster dataset
	 * @throws std::runtime_error if unable to get geotransform
	 */
	void createFromDataset(GDALDataset *p_dataset);

	/**
	 * Deconstructor for GDALRasterWrapper class. This method calls
	 * CPLFree() on any allocated raster buffers.
	 */
	~GDALRasterWrapper();

	/**
	 * Getter method for wrapped dataset.
	 *
	 * @returns GDALDataset *pointer to the underlying dataset
	 */
	GDALDataset *getDataset();

	/**
	 * Getter method for the raster driver.
	 *
	 * @returns std::string of short and long names of the raster driver
	 */
	std::string getDriver();

	/**
	 * Getter method for the coordinate reference system.
	 *
	 * @returns std::string json string of the CRS.
	 * @throws std::runtime_error if unable to acquire CRS
	 */
	std::string getCRS();

	/**
	 * Getter method for the raster width.
	 *
	 * @returns int raster width (x)
	 */
	int getWidth();

	/**
	 * Getter method for the raster height.
	 *
	 * @returns int raster height (y)
	 */
	int getHeight();

	/**
	 * Getter method for the number of raster bands.
	 *
	 * @returns int number of raster bands
	 */
	int getBandCount();

	/**
	 * Getter method for the maximum x value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double max x value
	 */
	double getXMax();

	/**
	 * Getter method for the minimum x value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double min x value
	 */
	double getXMin();

	/**
	 * Getter method for the maximum y value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double max y value
	 */
	double getYMax();

	/**
	 * Getter method for the minimum y value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double min y value
	 */
	double getYMin();

	/**
	 * Getter method for the pixel width. Scalar (absolute) value is given.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double pixel width
	 */
	double getPixelWidth();

	/**
	 * Getter method for the pixel height. Scalar (absolute) value is given.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double pixel height
	 */
	double getPixelHeight();

	/**
	 * Getter method for raster band names. Bands occur in order, meaning
	 * bands[0] corresponds to band 1, bands[1] to band 2, etc.
	 *
	 * @returns std::vector<std::string> vector of raster band names
	 */
	std::vector<std::string> getBands();	

	/**
	 * Getter method for geotransform.
	 *
	 * @returns double *array of 6 doubles representing GDAL geotransform
	 */
	double *getGeotransform();

	/**
	 * Getter method for a specific bands nodata value.
	 *
	 * @returs double nodata value
	 */
	double getBandNoDataValue(int band);

	/**
	 * Getter method for the raster image, used by the Python side 
	 * of the application. This function allocates and reads a raster band if necessary, 
	 * and uses py::memoryview::from_buffer() to create the buffer of the 
	 * correct size/dimensions without copying data unecessarily.
	 *
	 * This function requires that width and height be defined according
	 * to GDAL target_downscaling_factor rules. Otherwise, the incorrect amount
	 * of memory will be allocated. Information on target_downsampling_factor
	 * can be found here:
	 * https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec
	 *
	 * for py::memoryview::from_buffer() information see:
	 * https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	 *
	 * @param int width
	 * @param int height
	 * @param int band
	 * @throws std::runtime_error if unable to read raster band during allocation
	 * @returns py::buffer python memoryview of the raster
	 */
	py::buffer getRasterBandAsMemView(int width, int height, int band);

	/**
	 * Getter method for a GDALRasterBand in the raster, used by the C++ side of the application.
	 *
	 * @param int band zero-indexed
	 * @returns void *pointer to band
	 * @throws std::runtime_error if unable to read raster band during allocation
	 */
	GDALRasterBand *getRasterBand(int band);

	/**
	 * Getter method for the pixel / raster data type.
	 *
	 * @param int band to get type from
	 * @returns GDALDataType the data type
	 */
	GDALDataType getRasterBandType(int band);

	/**
	 * Getter method for the pixel /raster data type size.
	 *
	 * @param int band to get type size from
	 * @returns size_t the data type size in bytes
	 */
	size_t getRasterBandTypeSize(int band);

	/**
	 * Writes the raster to a specific file given by filename, by creating
	 * a copy of the GDALDatset.
	 *
	 * @param std::string filename
	 */
	void write(std::string filename);
};
