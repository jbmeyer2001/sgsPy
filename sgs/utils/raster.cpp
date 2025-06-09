/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for raster operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include "raster.h"

/******************************************************************************
			      GDALRasterWrapper()
******************************************************************************/
GDALRasterWrapper::GDALRasterWrapper(std::string filename) {
	//must register drivers before trying to open a dataset
	GDALAllRegister();

	//dataset
	this->p_dataset = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly)));
	if (!this->p_dataset) {
		throw std::runtime_error("dataset pointer is null after initialization, dataset unable to be initialized.");
	}

	//geotransform
	CPLErr cplerr = this->p_dataset->GetGeoTransform(this->geotransform);
	if (cplerr) {
		throw std::runtime_error("error getting geotransform from dataset.");
	}	
}

/******************************************************************************
			      ~GDALRasterWrapper()
******************************************************************************/
GDALRasterWrapper::~GDALRasterWrapper() {
	if (this->p_fullRaster) {
		CPLFree(this->p_fullRaster);
	}

	if (this->p_downsampledRaster) {
		CPLFree(this->p_downsampledRaster);
	}
}

/******************************************************************************
				  getDriver()
******************************************************************************/
std::string GDALRasterWrapper::getDriver() {
	return std::string(this->p_dataset->GetDriverName()) 
		+ "/" 
		+ std::string(this->p_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
}

/******************************************************************************
				    getCRS()
******************************************************************************/
std::string GDALRasterWrapper::getCRS() {
	char *p_crs;

	OGRErr ogrerr = OGRSpatialReference(this->p_dataset->GetProjectionRef()).exportToPROJJSON(&p_crs, nullptr);
	if (ogrerr) {
		throw std::runtime_error("error getting coordinate reference system from dataset.");
	}
		
	return std::string(p_crs);
}

/******************************************************************************
				   getWidth()
******************************************************************************/
int GDALRasterWrapper::getWidth() {
	return this->p_dataset->GetRasterXSize();
}

/******************************************************************************
				  getHeight()
******************************************************************************/
int GDALRasterWrapper::getHeight() {
	return this->p_dataset->GetRasterYSize();
}

/*****************************************************************************
				 getBandCount()
******************************************************************************/
int GDALRasterWrapper::getBandCount() {
	return this->p_dataset->GetRasterCount();
}

/******************************************************************************
				   getXMax()
******************************************************************************/
double GDALRasterWrapper::getXMax() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::max(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
	);
}

/******************************************************************************
				   getXMin()
******************************************************************************/
double GDALRasterWrapper::getXMin() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::min(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
	);
}

/******************************************************************************
				   getYMax()
******************************************************************************/
double GDALRasterWrapper::getYMax() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::max(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
	);
}

/******************************************************************************
				   getYMin()
******************************************************************************/
double GDALRasterWrapper::getYMin() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::min(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
	);
}

std::vector<double> GDALRasterWrapper::getOrigin() {
	std::vector<double> origin = {this->geotransform[0], this->geotransform[3]};
	return origin;
}

/******************************************************************************
				getPixelWidth()
******************************************************************************/
double GDALRasterWrapper::getPixelWidth() {
	return std::abs(this->geotransform[5]);
}

/******************************************************************************
				getPixelHeight()
******************************************************************************/
double GDALRasterWrapper::getPixelHeight() {
	return std::abs(this->geotransform[1]);
}

/******************************************************************************
				   getBands()
******************************************************************************/
std::vector<std::string> GDALRasterWrapper::getBands() {
	std::vector<std::string> retval;

	for( auto&& p_band : this->p_dataset->GetBands() ){
		retval.push_back(p_band->GetDescription());
	}

	return retval;
}

/******************************************************************************
				allocateRaster()
******************************************************************************/
void *GDALRasterWrapper::allocateRaster(int width, int height) {
	//get type and size information
	GDALDataType type = this->p_dataset->GetRasterBand(1)->GetRasterDataType();
	size_t size = GDALExtendedDataType::Create(type).GetSize();

	//allocate the whole raster
	void *p_raster = CPLMalloc(size * this->getBandCount() * width * height);

	//calculate the number of bytes that a single raster layer/band takes
	size_t layerSize = size * width * height;

	//iterate and read the raster bands. Begin at the start of the allocated chunk.
	void *bufferLocation = p_raster;
	for( auto&& p_band : this->p_dataset->GetBands() ) {
		/**
		 * NOTE: RasterIO should take care of windowing for us, so we can allocate the
		 * whole raster, with all bands, and let GDAL worry about RAM size, etc.
		 * 
		 * see https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec
		 * and https://gdal.org/en/stable/tutorials/raster_api_tut.html#reading-raster-data
		 */

		//perform raster read on current band
		CPLErr err = p_band->RasterIO(
			GF_Read, 			//GDALRWFlag eRWFlag
			0, 					//int nXOff
			0,					//int nYOff
			this->getWidth(),	//int nXSize
			this->getHeight(),	//int nYSize
			bufferLocation,		//void *pData
			width,				//int nBufXSize
			height,				//int nBufYSize
			type,				//GDALDataType eBufType 
			0,					//int nPixelSpace
			0					//int nLineSpace
		);
		if (err) {
			throw std::runtime_error("error reading raster band from dataset.");
		}

		//move the write buffer by the size of a band
		bufferLocation = (void *)((size_t)bufferLocation + layerSize);
	}

	return p_raster;
}

/******************************************************************************
				  getBuffer()
******************************************************************************/
template <typename T> 
py::buffer GDALRasterWrapper::getBuffer(size_t size, void *p_raster, int width, int height) {
	//see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	return py::memoryview::from_buffer(
		(T*)p_raster, 								//buffer
		{this->getBandCount(), height, width}, 		//shape
		{size * height * width, size * width, size} //stride
	);
}	

/******************************************************************************
			      getRasterAsMemView()
******************************************************************************/
py::buffer GDALRasterWrapper::getRasterAsMemView(int width, int height) {
	bool downsampled = (width != this->getWidth() || height != this->getHeight());

	//allocate full raster if required
	if (!downsampled && !this->fullRasterAllocated) {
		this->p_fullRaster = this->allocateRaster(width, height);
		this->fullRasterAllocated = true;
	}
	
	//allocate (or re-allocate) downsampled raster if required
	if (downsampled && (width != this->downsampledRasterWidth || height != this->downsampledRasterHeight)) {
		if (this->downsampledRasterAllocated) {
			CPLFree(this->p_downsampledRaster);
		}
		this->p_downsampledRaster = this->allocateRaster(width, height);
		this->downsampledRasterAllocated = true;
		this->downsampledRasterHeight = height;
		this->downsampledRasterWidth = width;
	}

	void *p_raster = downsampled ? this->p_downsampledRaster : this->p_fullRaster;	
	
	switch(this->p_dataset->GetRasterBand(1)->GetRasterDataType()) {
		case GDT_Int8:
			return getBuffer<int8_t>(1, p_raster, width, height);
		case GDT_UInt16:
			return getBuffer<uint16_t>(2, p_raster, width, height);
		case GDT_Int16:
			return getBuffer<int16_t>(2, p_raster, width, height);
		case GDT_UInt32:
			return getBuffer<uint32_t>(4, p_raster, width, height);
		case GDT_Int32:
			return getBuffer<int32_t>(4, p_raster, width, height);
		case GDT_Float32:
			return getBuffer<float>(4, p_raster, width, height);
		case GDT_UInt64:
			return getBuffer<uint64_t>(8, p_raster, width, height);
		case GDT_Int64:
			return getBuffer<int64_t>(8, p_raster, width, height);
		case GDT_Float64:
			return getBuffer<double>(8, p_raster, width, height);
		default:
			throw std::runtime_error("raster pixel data type not supported.");
	}
}

/******************************************************************************
				  getRaster()				     
******************************************************************************/
void *GDALRasterWrapper::getRaster() {
	if (!this->fullRasterAllocated) {
		this->p_fullRaster = this->allocateRaster(this->getWidth(), this->getHeight());
		this->fullRasterAllocated = true;
	}

	return this->p_fullRaster;
}
