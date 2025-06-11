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

	//noData representation
	switch (this->p_dataset->GetRasterBand(1)->GetRasterDataType()) {
		case GDT_Int64:
			this->noDataValue<int64_t> = this->p_dataset->GetRasterBand(1)->GetNoDataValueAsInt64();
			break;
		case GDT_UInt64:
			this->noDataValue<uint64_t> = this->p_dataset->GetRasterBand(1)->GetNoDataValueAsUInt64();
			break;
		default:
			this->noDataValue<double> = this->p_dataset->GetRasterBand(1)->GetNoDataValue();
	}

	switch(this->p_dataset->GetRasterBand(1)->GetRasterDataType()) {
		case GDT_Int8:
			this->pixelSize = 1;
			break;
		case GDT_UInt16:
		case GDT_Int16:
			this->pixelSize = 2;
			break;
		case GDT_UInt32:
		case GDT_Int32:
		case GDT_Float32:
			this->pixelSize = 4;
			break;
		case GDT_UInt64:
		case GDT_Int64:
		case GDT_Float64:
			this->pixelSize = 8;
			break;
		default:
			throw std::runtime_error("raster pixel data type not supported.");
	}

	//TESTING CODE FOR NODATA ADJUSTMENT FUNCTIONALITY
	//
	//THIS WILL BE DELETED
	//
	//IN THE FUTURE, A MORE CONCRETE TESTING STRATEGY WILL HAVE TO BE EMPLOYED TO ENSURE CHANGES
	//TO THE GDALRASTERWRAPPER CLASS DONT IMPEDE THE NODATA ADJUSTMENT FUNCTIONALITY
	
	std::cout << "RASTER:" << std::endl;
	double *raster = (double *)this->getRaster();
	for (size_t i = 0; i < this->getWidth() * this->getHeight(); i++) {
		std::cout << "[" << i << "] = " << raster[i] << std::endl;
	}

	std::cout << std::endl << std::endl << std::endl << "ADJUSTED RASTER: " << std::endl;
	raster = (double *)this->getNoDataRaster();
	for (size_t i = 0; i < this->getWidth() * this->getHeight() - this->noDataCount; i++) {
		std::cout << "[" << i << "] = " << raster[i] << std::endl;
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
size_t GDALRasterWrapper::getWidth() {
	return (size_t)this->p_dataset->GetRasterXSize();
}

/******************************************************************************
				  getHeight()
******************************************************************************/
size_t GDALRasterWrapper::getHeight() {
	return (size_t)this->p_dataset->GetRasterYSize();
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
	return std::max(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * this->getWidth() + this->geotransform[2] * this->getHeight()
	);
}

/******************************************************************************
				   getXMin()
******************************************************************************/
double GDALRasterWrapper::getXMin() {
	return std::min(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * this->getWidth() + this->geotransform[2] * this->getHeight()
	);
}

/******************************************************************************
				   getYMax()
******************************************************************************/
double GDALRasterWrapper::getYMax() {
	return std::max(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * this->getWidth() + this->geotransform[5] * this->getHeight()
	);
}

/******************************************************************************
				   getYMin()
******************************************************************************/
double GDALRasterWrapper::getYMin() {
	return std::min(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * this->getWidth() + this->geotransform[5] * this->getHeight()
	);
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
void *GDALRasterWrapper::allocateRaster(size_t width, size_t height) {
	//get type and size information
	GDALDataType type = this->p_dataset->GetRasterBand(1)->GetRasterDataType();
	size_t size = GDALExtendedDataType::Create(type).GetSize();

	//allocate the whole raster
	void *p_raster = CPLMalloc(size * this->getBandCount() * width * height);

	//calculate the number of bytes that a single raster layer/band takes
	size_t layerSize = (size_t)size * width * height;

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
py::buffer GDALRasterWrapper::getBuffer(size_t size, void *p_raster, size_t width, size_t height) {
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
py::buffer GDALRasterWrapper::getRasterAsMemView(size_t width, size_t height) {
	//allocate (or re-allocate) display raster if required
	if (width != this->displayRasterWidth || height != this->displayRasterHeight) {
		if (this->downsampledRasterAllocated) {
			CPLFree(this->p_downsampledRaster);
		}
		this->p_downsampledRaster = this->allocateRaster(width, height);
		this->displayRasterAllocated = true;
		this->displayRasterHeight = height;
		this->displayRasterWidth = width;
	}
	
	switch(this->p_dataset->GetRasterBand(1)->GetRasterDataType()) {
		case GDT_Int8:
			return getBuffer<int8_t>(1, this->p_displayRaster, width, height);
		case GDT_UInt16:
			return getBuffer<uint16_t>(2, this->p_displayRaster, width, height);
		case GDT_Int16:
			return getBuffer<int16_t>(2, this->p_displayRaster, width, height);
		case GDT_UInt32:
			return getBuffer<uint32_t>(4, this->p_displayRaster, width, height);
		case GDT_Int32:
			return getBuffer<int32_t>(4, this->p_displayRaster, width, height);
		case GDT_Float32:
			return getBuffer<float>(4, this->p_displayRaster, width, height);
		case GDT_UInt64:
			return getBuffer<uint64_t>(8, this->p_displayRaster, width, height);
		case GDT_Int64:
			return getBuffer<int64_t>(8,  this->p_displayRaster, width, height);
		case GDT_Float64:
			return getBuffer<double>(8, this->p_displayRaster, width, height);
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

/******************************************************************************
				   isNoData()				     
******************************************************************************/
inline bool isNoData(double val) { return std::is_nan(val) || val == this->noDataValue; }
inline bool isNoData(int64_t val) { return std::is_nan(val) || val == this->noDataValue; }
inline bool isNoData(uint64_t val) { return std::is_nan(val) || val == this->noDataValue; }


/******************************************************************************
				 copyInChunks()				     
******************************************************************************/
void GDALRasterWrapper::copyInChunks(void *memcpySrc, void *memcpyDst, size_t memcpySize) {
	size_t memcpyGapSize = (size_t)(dataBlockStart - noDataBlockStart);

	size_t bytesCopied = 0;
	while (bytesCopied < memcpySize) {
		size_t copyAmount = std::min(memcpyGapSize, memcpySize - bytesCopied);
		std::memcpy(
			(void *)((size_t)memcpyDst + bytesCopied),	//destination
			(void *)((size_t)memcpySrc + bytesCopied),	//source
			copyAmount									//bytes to copy
		);
		bytesCopied += copyAmount;
	}
}

/******************************************************************************
			       getNoDataRaster()				     
******************************************************************************/
void *GDALRasterWrapper::getNoDataRaster() {
	if (this->nanAdjusted) {
		return p_raster;
	}

	this->getRaster();

	void *noDataBlockStart = nullptr;
	void *dataBlockStart = nullptr;
	bool prevNoData = isNoData(p_raster[0]); 
	bool curNoData;

	if (prevNoData) {
		noDataBlockStart = this->p_raster;
		this->noDataCount++;
	}
	else {
		dataBlockStart = this->p_raster;
		originalIndex.push_back(0);
		adjustedIndex.push_back(0);
	}
	
	for (size_t i = 1; i < this->getWidth() * this->getHeight(); i++) {
		bool curNoData = isNoData(this->p_raster[i]);
		
		if(curNoData && prevNoData) {
			this->noDataCount++;
		}
		else if (!curNoData && prevNoData) {
			dataBlockStart = (void *)((size_t)this->p_raster + (i * this->pixelSize))
			originalIndex.push_back(i);
			adjustedIndex.push_back(i - this->noDataCount);
		}
		else if (curNoData && !prevNoData) {
			this->noDataCount++;			
			
			if (!noDataBlockStart) {
				noDataBlockStart = (void *)((size_t)this->p_raster + (i * this->pixelSize));
				prevNoData = curNoData;
				continue;
			}

			size_t dataBlockSize = (size_t)this->p_raster + ((i - 1) * this->pixelSize) - (size_t)dataBlockStart;
			copyInChunks(dataBlockStart, noDataBlockStart, dataBlockSize)

			noDataBlockStart = (void *)((size_t)noDataBlockStart + memcpySize);
			dataBlockStart = nullptr;
		}
		
		prevNoData = curNoData;
	}

	//if the last pixel was a data pixel, and there were any noData pixels in the image, we have
	//one more block to copy
	if (this->noDataCount > 0 && !prevNoData) {
		size_t dataBlockSize = (siz_t)this->p_raster + ((this->getHeight() * this->getWidth() - 1) * this->size) - (size_t)dataBlockStart;
		copyInChunks(dataBlockStart, noDataBlockStart, dataBlockSize)
	}
}

/******************************************************************************
			       getOriginalIndex()				     
******************************************************************************/
size_t GDALRasterWrapper::getOriginalIndex(size_t adjustedIndex) {
	//apparently according to cppreference:
	//std::upper_bound(x) returns a value such that x < val
	//std::lower_bound(x) returns a value such that x <= val
	//
	//Make it  make sense??? Shouldn't lower_bound return x > val ???
	//
	//Anyways... 
	//We want to get the first value smaller than or equal to x. This can be done by using upper_bound(x) - 1.
	//
	//Using that we can get the number of noData pixels which occured before the adjusted index (the adjustment), which
	//when added to the adjusted index gives the original index.
	auto boundIt = std::upper_bound(this->adjustedIndexes.begin(), this->adjustedIndexes.end(), adjustedIndex);
	
	if (boundIt == this->adjustedIndexes.end()) {
		size_t boundIndex = this->adjustedIndexes.back();
	}
	else {
		size_t boundIndex = std::distance(this->adjustedIndexes.begin(), boundIt) - 1;
	}

	size_t adjustemnt = this->originalIndexes.at(boundIndex) - this->adjustedIndexes(boundIndex);
	return adjustedIndex + adjustment;
}
