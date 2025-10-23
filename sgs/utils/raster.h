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
	std::string crs = "";
	char *p_proj = nullptr;

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
	void readRasterBand(int width, int height, int band) {	
		GDALDataType type = this->getRasterBandType(band);
		size_t size = this->getRasterBandTypeSize(band);
		size_t max = std::numeric_limits<size_t>::max();
		
		//perform size checks
		if (max / size < static_cast<size_t>(width) ||
		    max / (size * static_cast<size_t>(width)) < static_cast<size_t>(height)) {
				throw std::runtime_error("raster too large to fit in memory.");
		}

		if (static_cast<size_t>(height) * static_cast<size_t>(width) * size > GIGABYTE) {
				throw std::runtime_error("sgs does not allow allocation of a raster into memory for direct pixel access purposes if it would be larger than 1 gigabyte.");
		}

		//allocate data
		void *p_data = VSIMalloc3(height, width, size);

		//perform raster read on current band
		CPLErr err = this->p_dataset->GetRasterBand(band + 1)->RasterIO(
			GF_Read, 			//GDALRWFlag eRWFlag
			0, 				//int nXOff
			0,				//int nYOff
			this->getWidth(),		//int nXSize
			this->getHeight(),		//int nYSize
			p_data,				//void *pData
			width,				//int nBufXSize
			height,				//int nBufYSize
			type,				//GDALDataType eBufType 
			0,				//int nPixelSpace
			0				//int nLineSpace
		);
		if (err) {
			throw std::runtime_error("error reading raster band from dataset.");
		}

		//update dislpay information as required
		if (width != this->getWidth() || height != this->getHeight()) {
			this->displayRasterBandRead[band] = true;
			this->displayRasterBandPointers[band] = p_data;
		}
		else {
			this->rasterBandRead[band] = true;
			this->rasterBandPointers[band] = p_data;
		}
	}

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
	py::buffer getBuffer(size_t size, void *p_buffer, int width, int height) {
		//see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
		return py::memoryview::from_buffer(
			(T*)p_buffer, 		//buffer
			{height, width}, 	//shape
			{size * width, size} 	//stride
		);
	}	
	
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
	GDALRasterWrapper(std::string filename) {
		//must register drivers before trying to open a dataset
		GDALAllRegister();

		//dataset
		GDALDataset *p_dataset = GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly));
		if (!p_dataset) {
			throw std::runtime_error("dataset pointer is null after initialization, dataset unable to be initialized.");
		}

		this->createFromDataset(p_dataset);
	}	

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
	GDALRasterWrapper(GDALDataset *p_dataset, std::vector<void *> bands) {
		this->createFromDataset(p_dataset);
		this->rasterBandPointers = bands;
		this->rasterBandRead = std::vector<bool>(bands.size(), true);
	}

	/**
	 * Constructor for GDALRasterWrapper class using just a GDAL 
	 * dataset pointer, by calling createFromDataset().
	 *
	 * @param GDALDataset *p_dataset GDAL raster dataset
	 */
	GDALRasterWrapper(GDALDataset *p_dataset) {
		this->createFromDataset(p_dataset);
	}

	/**
	 * Populates/constructs the GDALRasterWrapper object using a raster
	 * dataset pointer. Called by some of the GDALRasterWrapper 
	 * constructors as a helper function.
	 *
	 * @param GDALDataset *GDAL raster dataset
	 * @throws std::runtime_error if unable to get geotransform
	 */
	void createFromDataset(GDALDataset *p_dataset) {
		this->p_dataset = GDALDatasetUniquePtr(p_dataset);

		//geotransform
		CPLErr cplerr = this->p_dataset->GetGeoTransform(this->geotransform);
		if (cplerr) {
			throw std::runtime_error("error getting geotransform from dataset.");
		}

		//crs
		this->crs = std::string(OGRSpatialReference(this->p_dataset->GetProjectionRef()).GetName());

		//initialize (but don't read) raster band pointers
		this->rasterBandPointers = std::vector<void *>(this->getBandCount(), nullptr);
		this->rasterBandRead = std::vector<bool>(this->getBandCount(), false);
		this->displayRasterBandPointers = std::vector<void *>(this->getBandCount(), nullptr);
		this->displayRasterBandRead = std::vector<bool>(this->getBandCount(), false);	
	}

	/**
	 * Deconstructor for GDALRasterWrapper class. This method calls
	 * CPLFree() on any allocated raster buffers.
	 */
	~GDALRasterWrapper() {
		for (int i = 0; i < this->getBandCount(); i++) {
			if (this->rasterBandRead[i]) {
				CPLFree(this->rasterBandPointers[i]);
			}

			if (this->displayRasterBandRead[i]) {
				CPLFree(this->displayRasterBandPointers[i]);
			}
		}

		if (this->p_proj) {
			free(this->p_proj);
		}	
	}

	/**
	 * Getter method for wrapped dataset.
	 *
	 * @returns GDALDataset *pointer to the underlying dataset
	 */
	GDALDataset *getDataset() {
		return this->p_dataset.get();
	}

	/**
	 * Getter method for the raster driver.
	 *
	 * @returns std::string of short and long names of the raster driver
	 */
	std::string getDriver() {
		return std::string(this->p_dataset->GetDriverName()) 
			+ "/" 
			+ std::string(this->p_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
	}

	/* Getter method for the full projection information as wkt.
	 *
	 * @param std::string projection as wkt
	 */
	std::string getFullProjectionInfo() {
		if (!this->p_proj) {
			OGRErr ogrerr = OGRSpatialReference(this->p_dataset->GetProjectionRef()).exportToPrettyWkt(&this->p_proj);
			if (ogrerr) {
				throw std::runtime_error("error getting projection as WKT from dataset.");
			}
		}

		return std::string(this->p_proj);
	}
	/**
	 * Get the CRS name from the OGRSpatialReference object
	 *
	 * @returns CRS name
	 */
	std::string getCRS(){
		return this->crs;
	}

	/**
	 * Getter method for the raster width.
	 *
	 * @returns int raster width (x)
	 */
	int getWidth() {
		return this->p_dataset->GetRasterXSize();
 	}

	/**
	 * Getter method for the raster height.
	 *
	 * @returns int raster height (y)
	 */
	int getHeight(){
		return this->p_dataset->GetRasterYSize();
	}

	/**
	 * Getter method for the number of raster bands.
	 *
	 * @returns int number of raster bands
	 */
	int getBandCount() {
		return this->p_dataset->GetRasterCount();
	}

	/**
	 * Getter method for the maximum x value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double max x value
	 */
	double getXMax() {
		int width = this->p_dataset->GetRasterXSize();
		int height = this->p_dataset->GetRasterYSize();
		return std::max(
			this->geotransform[0], 
			this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
		);
	}

	/**
	 * Getter method for the minimum x value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double min x value
	 */
	double getXMin() {
		int width = this->p_dataset->GetRasterXSize();
		int height = this->p_dataset->GetRasterYSize();
		return std::min(
			this->geotransform[0], 
			this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
		);
	}

	/**
	 * Getter method for the maximum y value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double max y value
	 */
	double getYMax() {
		int width = this->p_dataset->GetRasterXSize();
		int height = this->p_dataset->GetRasterYSize();
		return std::max(
			this->geotransform[3],
			this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
		);
	}

	/**
	 * Getter method for the minimum y value in georeferenced coordinate space.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double min y value
	 */
	double getYMin() {
		int width = this->p_dataset->GetRasterXSize();
		int height = this->p_dataset->GetRasterYSize();
		return std::min(
			this->geotransform[3],
			this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
		);
	}

	/**
	 * Getter method for the pixel width. Scalar (absolute) value is given.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double pixel width
	 */
	double getPixelWidth() {
		return std::abs(this->geotransform[5]);
	}

	/**
	 * Getter method for the pixel height. Scalar (absolute) value is given.
	 * see https://gdal.org/en/stable/tutorials/geotransforms_tut.html
	 *
	 * @returns double pixel height
	 */
	double getPixelHeight() {
		return std::abs(this->geotransform[1]);
	}

	/**
	 * Getter method for raster band names. Bands occur in order, meaning
	 * bands[0] corresponds to band 1, bands[1] to band 2, etc.
	 *
	 * @returns std::vector<std::string> vector of raster band names
	 */
	std::vector<std::string> getBands(){
		std::vector<std::string> retval;

		for( auto&& p_band : this->p_dataset->GetBands() ){
			retval.push_back(p_band->GetDescription());
		}

		return retval;
	}	

	/**
	 * Getter method for geotransform.
	 *
	 * @returns double *array of 6 doubles representing GDAL geotransform
	 */
	double *getGeotransform() {
		return this->geotransform;
	}

	/**
	 * Getter method for a specific (0-indexed) bands nodata value.
	 *
	 * @returs double nodata value
	 */
	double getBandNoDataValue(int band){
		GDALRasterBand *p_band = this->getRasterBand(band);
		return p_band->GetNoDataValue();
	}

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
	py::buffer getRasterBandAsMemView(int width, int height, int band) {
		bool display = (width != this->getWidth() || height != this->getHeight());
		void *p_buffer;
		GDALDataType type = this->getRasterBandType(band);

		//allocate raster if required
		if (!display && !this->rasterBandRead[band]) {
			this->readRasterBand(width, height, band);
		}

		//(re)allocate display raster if required
		if (display) {
			if (width != this->displayRasterWidth || height != this->displayRasterHeight) {
				free(this->displayRasterBandPointers[band]);
				this->displayRasterBandRead[band] = false;
			}

			if (this->displayRasterBandRead[band] == false) {
				this->readRasterBand(width, height, band);
			}
		}

		//get the (allocated) data buffer
		p_buffer = (!display) ?
			this->rasterBandPointers[band] :
			this->displayRasterBandPointers[band];

		switch(type) {
			case GDT_Int8:
				return getBuffer<int8_t>(sizeof(int8_t), p_buffer, width, height);
			case GDT_UInt16:
				return getBuffer<uint16_t>(sizeof(uint16_t), p_buffer, width, height);
			case GDT_Int16:
				return getBuffer<int16_t>(sizeof(int16_t), p_buffer, width, height);
			case GDT_UInt32:
				return getBuffer<uint32_t>(sizeof(uint32_t), p_buffer, width, height);
			case GDT_Int32:
				return getBuffer<int32_t>(sizeof(int32_t), p_buffer, width, height);
			case GDT_Float32:
				return getBuffer<float>(sizeof(float), p_buffer, width, height);
			case GDT_Float64:
				return getBuffer<double>(sizeof(double), p_buffer, width, height);
			default:
				throw std::runtime_error("raster pixel data type not supported.");
		}
	}

	/**
	 * Getter method for a GDALRasterBand in the raster, used by the C++ side of the application.
	 *
	 * @param int band zero-indexed
	 * @returns GDALRasterBand *pointer to band
	 * @throws std::runtime_error if unable to read raster band during allocation
	 */
	GDALRasterBand *getRasterBand(int band) {
		return this->p_dataset->GetRasterBand(band + 1);	
	}

	/**
	 * Getter method for the whole GDALRasterBand data buffer.
	 *
	 * @param int band zero-indexed
	 * @returns void *data pointer to band
	 */
	void *getRasterBandBuffer(int band) {
		if (!this->rasterBandRead[band]) {
			this->readRasterBand(
				this->getWidth(), 
				this->getHeight(),
				band
			);
		}

		return this->rasterBandPointers[band];
	}

	/**
	 * Getter method for the pixel / raster data type.
	 *
	 * @param int band to get type from
	 * @returns GDALDataType the data type
	 */
	GDALDataType getRasterBandType(int band) {
		GDALRasterBand *p_band = this->p_dataset->GetRasterBand(band + 1);
		return p_band->GetRasterDataType();	
	}

	/**
	 * Getter method for the pixel /raster data type size.
	 *
	 * @param int band to get type size from
	 * @returns size_t the data type size in bytes
	 */
	size_t getRasterBandTypeSize(int band) {
		switch (this->getRasterBandType(band)) {
			case GDALDataType::GDT_Int8:
				return 1;
			case GDALDataType::GDT_UInt16:
			case GDALDataType::GDT_Int16:
				return 2;
			case GDALDataType::GDT_UInt32:
			case GDALDataType::GDT_Int32:
			case GDALDataType::GDT_Float32:
				return 4;
			case GDALDataType::GDT_Float64:
				return 8;
			default:
				std::string errorMsg = "GDALDataType of band " + std::to_string(band) + " not supported.";
				throw std::runtime_error(errorMsg);
		}
	}

	/**
	 * Writes the raster to a specific file given by filename, by creating
	 * a copy of the GDALDatset.
	 *
	 * @param std::string filename
	 */
	void write(std::string filename) {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();
		
		if (extension != ".tif") {
			throw std::runtime_error("write only supports .tif files right now");
		}
	
		GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	
		GDALClose(p_driver->CreateCopy(filename.c_str(), this->p_dataset.get(), (int)false, nullptr, nullptr, nullptr));
	}
};
