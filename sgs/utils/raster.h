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
	void *p_raster = nullptr;
	py::buffer numpyRaster;
	bool rasterAllocated = false;
	double geotransform[6];

	/**
	 * Internal function used to allocate the raster data. Memory
	 * is allocated using CPLMalloc(). Rasters are read using
	 * GDALDataset::RasterIO().
	 *
	 * @throws std::runtime_error if unable to read raster band
	 */
	void allocateRaster();

	/**
	 * Internal function which returns a pybuffer of the raster, using
	 * the type specified.
	 */
	template <typename T> py::buffer getBuffer(size_t size);
	
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
	 * CPLFree() on the raster pointer if memory was allocated.
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
	 * Getter method for the raster image, used by the Python side 
	 * of the application. This function uses py::memoryview::from_buffer() 
	 * to create the buffer of the correct size/dimensions.
	 *
	 * The memory view is then used to initialize a NumPy array without
	 * copying the data.
	 *
	 * see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	 *
	 * @throws std::runtime_error if unable to read raster band during allocation
	 */
	py::buffer getRasterAsMemView();

	/**
	 * Getter method for the raster image, used by the C++ side of the application.
	 *
	 * @throws std::runtime_error if unable to read raster band during allocation
	 */
	void *getRaster();
};

/** 
 * pybind11 module for exposing GDALRasterWrapper data to Python.
 * Define only constructor and relevant getter functions.
 */
PYBIND11_MODULE(raster, m) {
	py::class_<GDALRasterWrapper>(m, "GDALRasterWrapper")
		.def(py::init<std::string>())
		.def("get_driver", &GDALRasterWrapper::getDriver)
		.def("get_crs", &GDALRasterWrapper::getCRS)
		.def("get_height", &GDALRasterWrapper::getHeight)
		.def("get_width", &GDALRasterWrapper::getWidth)
		.def("get_band_count", &GDALRasterWrapper::getBandCount)
		.def("get_xmin", &GDALRasterWrapper::getXMin)
		.def("get_xmax", &GDALRasterWrapper::getXMax)
		.def("get_ymin", &GDALRasterWrapper::getYMin)
		.def("get_ymax", &GDALRasterWrapper::getYMax)
		.def("get_pixel_height", &GDALRasterWrapper::getPixelHeight)
		.def("get_pixel_width", &GDALRasterWrapper::getPixelWidth)
		.def("get_bands", &GDALRasterWrapper::getBands)
		.def("get_raster_as_memoryview", &GDALRasterWrapper::getRasterAsMemView);
}
