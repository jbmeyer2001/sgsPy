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
 * In Linux, memory will be virtually allocated using GDALs existing functionality 
 * to allow for large images to be used without tiling or the use of scanlines. 
 * Windows functionality has not been decided, as GDAL virtual memory only exists 
 * on linux currently.
 */
class GDALRasterWrapper {
	private:
	GDALDatasetUniquePtr p_dataset;
	void *p_raster;
	py::buffer numpyRaster;
	bool rasterAllocated = false;
	double geotransform[6];
	GDALDataType type;

	/**
	 * TODO document
	 */
	void allocateRaster();

	/**
	 * TODO document
	 */
	void allocateRasterHelper(size_t size);

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
	 * @param filename an std::string
	 */
	GDALRasterWrapper(std::string filename);

	/**
	 * Deconstructor for GDALRasterWrapper class. This method calls
	 * CPLVirtualMemFree() if virtual memory has been allocated.
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
	 * @returns python dict of the CRS.
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
	 * Getter method for the number of raster bands/layers.
	 *
	 * @returns int number of raster layers
	 */
	int getLayers();

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
	 * bands[0] corrosponds to band 1, bands[1] to band 2, etc.
	 *
	 * @returns std::vector<std::string> vector of raster band names.
	 */
	std::vector<std::string> getBands();	

	/**
	 * Getter method for the raster image in virtual memory, used
	 * by the python side of the application. This function
	 * uses the getBuffer() function to define the type of the
	 * buffer.
	 *
	 * Python memory view is used to create a numpy array.
	 *
	 * see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	 */
	py::buffer getRasterAsMemView();

	/**
	 * Getter method for the raster image in virtual memory, used
	 * by the C++ side of the application.
	 *
	 * This is a public method which used to expose the getRasterPointer()
	 * functionality to external C++ functions.
	 */
	void *getRaster();
};

PYBIND11_MODULE(raster, m) {
	py::class_<GDALRasterWrapper>(m, "GDALRasterWrapper")
		.def(py::init<std::string>())
		.def("get_driver", &GDALRasterWrapper::getDriver)
		.def("get_crs", &GDALRasterWrapper::getCRS)
		.def("get_height", &GDALRasterWrapper::getHeight)
		.def("get_width", &GDALRasterWrapper::getWidth)
		.def("get_layers", &GDALRasterWrapper::getLayers)
		.def("get_xmin", &GDALRasterWrapper::getXMin)
		.def("get_xmax", &GDALRasterWrapper::getXMax)
		.def("get_ymin", &GDALRasterWrapper::getYMin)
		.def("get_ymax", &GDALRasterWrapper::getYMax)
		.def("get_pixel_height", &GDALRasterWrapper::getPixelHeight)
		.def("get_pixel_width", &GDALRasterWrapper::getPixelWidth)
		.def("get_bands", &GDALRasterWrapper::getBands)
		.def("get_raster_as_memoryview", &GDALRasterWrapper::getRasterAsMemView);
}
