/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for raster operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#pragma once

#include <gdal_priv.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

	size_t pixelSize;

	void *p_raster = nullptr;
	bool rasterAllocated = false;
	bool noDataAdjusted = false;
	size_t noDataCount = 0;
	template<typename T> noDataValue;
	std::vector<size_t> originalIndexes;
	std::vector<size_t> adjustedIndexes;

	void *p_displayRaster = nullptr;
	bool displayRasterAllocated = false;
	size_t displayRasterWidth = 0;
	size_t displayRasterHeight = 0;

	double geotransform[6];

	/**
	 * Internal function used to allocate the raster data.
	 * This function is used for allocating both full and 
	 * downsampled rasters, which is why it takes the width and
	 * height parameters. 
	 *
	 * Memory is allocated using CPLMalloc(). Rasters are read 
	 * using GDALDataset::RasterIO().
	 *
	 * @param int width
	 * @param int height
	 * @returns void * allocated raster buffer
	 * @throws std::runtime_error if unable to read raster band
	 */
	void *allocateRaster(size_t width, size_t height);

	/**
	 * Internal function which returns a pybuffer of the display raster, 
	 * using the type specified. The raster must have already been allocated
	 * using allocateRaster(). 
	 *
	 * This function is used for getting the buffer to the display
	 * raster, which is why it takes the raster pointer,
	 * width, and height parameters. 
	 *
	 * @param size_t size of data type for one pixel
	 * @param void *p_raster pointer to allocated raster
	 * @param int width 
	 * @param int height
	 * @returns py::buffer memoryview object of the data
	 */
	template <typename T> 
	py::buffer getBuffer(size_t size, void *p_raster, size_t width, size_t height);

	/**
	 * Internal functions to check whether a given pixel value is no data.
	 *
	 * For GDT_Int64 and GDT_UInt64 images, using the typical GetNoDataValue()
	 * function which returns a double may be lossy, so int64_t and uint64_t
	 * values should be used. That is the purpose of the function overloads.
	 *
	 * @param double/int64_t/uint64_t val pixel vaue
	 * @returns true if noData, false if not noData
	 */
	inline bool isNoData(double val);
	inline bool isNoData(int64_t val);
	inline bool isNoData(uint64_t val);

	/**
	 * Internal helper function for getNoDataRaster(). 
	 *
	 * Due to the way
	 * that getNoDataRaster() works by copying the data sections to
	 * earlier indexes in the same raster previously occupied by noData,
	 * it's possible that there will be overlaps between the destination
	 * block and source block.
	 *
	 * When there is overlap between a src and dst chunk with std::memcpy
	 * there is undefined behavior. 
	 *
	 * This function copies the src data block to the index at memcpyDst
	 * in chunks to ensure there will be no overlap of src and dst in
	 * the call to std::memcpy.
	 */	
	void copyInChunks(void *memcpySrc, void *memcpyDst, size_t memcpySize);

	/**
	 * Internal function which returns a pointer to a noData adjusted version
	 * of the raster. The noData adjustd raster is the same raster except with
	 * no noData pixels. The length of the new raster will be: 
	 *
	 * height * width - num_noDataPixels
	 *
	 * This behavior will be required by many of the algorithms.
	 * TODO: explain more how function works, and how we can get the original index using the std::Vector
	 */
	void *getNoDataRaster();

	/**
	 * Takes the index of a pixel in the adjusted raster, and using the originalIndexes
	 * and adjustedIndexes vector, calculates what the index would have been in the original
	 * image before being noData adjusted.
	 *
	 * The function is necessary for getting the geographic coordinates of a particular adjusted
	 * pixel, since the geotransform calculations are based on the original image containing noData values.
	 *
	 * @param size_t the index of a pixel in the noData adjusted raster
	 * @param size_t the index that pixel would have been at in the original raster before noData adjustment.
	 */
	size_t getOriginalIndex(size_t adjustedIndex);

	public:
	/**
	 * Constructor for GDALRasterWrapper class. This method registers
	 * drivers, creates a GDALDataset object, and gets the geotransform
	 * information from the GDALDataset object.
	 *
	 * @param filename as std::string
	 * @throws std::runtime_error if dataset is not initialized
	 * @throws std::runtime_error if unable to get geotransform
	 */
	GDALRasterWrapper(std::string filename);

	/**
	 * Deconstructor for GDALRasterWrapper class. This method calls
	 * CPLFree() on any allocated raster buffers.
	 */
	~GDALRasterWrapper();

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
	 * Return is type size_t so when multiplying
	 * with height, there is no danger of accidental overflow.
	 *
	 * @returns size_t raster width (x)
	 */
	size_t getWidth();

	/**
	 * Getter method for the raster height.
	 * Return is type size_t so when multiplying
	 * with width, there is no danger of accidental overflow.
	 *
	 * @returns size_t raster height (y)
	 */
	size_t getHeight();

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
	 * Getter method for the raster image, used by the Python side 
	 * of the application. This function allocates the raster if necessary, and
	 * uses py::memoryview::from_buffer() to create the buffer of the 
	 * correct size/dimensions without copying data unecessarily.
	 *
	 * This function requires that width and height already be allocated according
	 * to GDAL target_downscaling_factor rules. Otherwise, the incorrect amount
	 * of memory will be allocated. Information on target_downsampling_factor
	 * can be found here:
	 * https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec
	 *
	 * A full raster will be allocated if the width and height correspond to the
	 * full raster size, and the p_fullRaster pointer has yet to be allocated.
	 * A downsampled raster will be allocated if a raster with the given width
	 * and height has yet to be allocated. Note, this may require de-allocation
	 * if a downsampled raster has been allocated with different width/height.
	 *
	 * for py::memoryview::from_buffer() information see:
	 * https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	 *
	 * @param int width
	 * @param int height
	 * @throws std::runtime_error if unable to read raster band during allocation
	 * @returns py::buffer python memoryview of the raster
	 */
	py::buffer getRasterAsMemView(int width, int height);

	/**
	 * Getter method for the raster image, used by the C++ side of the application.
	 *
	 * @throws std::runtime_error if unable to read raster band during allocation
	 */
	void *getRaster();
};
