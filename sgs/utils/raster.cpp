/******************************************************************************
 *
 * Project: sgs
 * Purpose: GDALDataset wrapper for raster operations
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include <iostream>

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

	//data type
	this->rasterType = GDALExtendedDataType::Create(this->p_dataset->GetRasterBand(1)->GetRasterDataType());

	//initialize (but don't read) raster band pointers
	this->rasterBandPointers = std::vector<void *>(this->getBandCount(), nullptr);
	this->rasterBandRead = std::vector<bool>(this->getBandCount(), false);
	this->displayRasterBandPointers = std::vector<void *>(this->getBandCount(), nullptr);
	this->displayRasterBandRead = std::vector<bool>(this->getBandCount(), false);	
}

/******************************************************************************
			      GDALRasterWrapper()
******************************************************************************/
GDALRasterWrapper::GDALRasterWrapper(
	std::vector<void *> rasterBands,
	std::vector<std::string> rasterBandNames,
	int width,
	int height,
	GDALDataType type,
	double *geotransform,
	std::string projection)
{
	//must register drivers before trying to open a dataset
	GDALAllRegister();

	//create in-memory dataset
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
	this->p_dataset = GDALDatasetUniquePtr(p_driver->Create(
		"",
		width,
		height,
		0,
		type,
		nullptr
	));
	if (!this->p_dataset) {
		throw std::runtime_error("dataset pointer is null after initialization, dataset unabel to be initialized.");
	}

	//set crs
	CPLErr err = this->p_dataset->SetProjection(projection.c_str());
	if (err) {
		throw std::runtime_error("error setting projection");
	}

	//dynamically add raster bands and their descriptions
	for (size_t i = 0; i < rasterBands.size(); i++) {
		char **papszOptions = nullptr;
		papszOptions = CSLSetNameValue(papszOptions, "DATAPOINTER", std::to_string((size_t)rasterBands[i]).c_str());
		err = this->p_dataset->AddBand(type, papszOptions);
		if (err) {
			throw std::runtime_error("error adding band.");
		}
		this->p_dataset->GetRasterBand(i + 1)->SetDescription(rasterBandNames[i].c_str());
	}	

	//copy geotransform values (not pointers) and set dataset geotransform
	this->geotransform[0] = geotransform[0];
	this->geotransform[1] = geotransform[1];
	this->geotransform[2] = geotransform[2];
	this->geotransform[3] = geotransform[3];
	this->geotransform[4] = geotransform[4];
	this->geotransform[5] = geotransform[5];
	err = this->p_dataset->SetGeoTransform(this->geotransform);
	if (err) {
		throw std::runtime_error("error setting geotransform.");
	}

	//set raster data type
	this->rasterType = GDALExtendedDataType::Create(type);

	//initialize raster band properties
	this->p_raster = rasterBandPointers[0];
	this->rasterAllocated = true;
	this->rasterBandPointers = rasterBands;
	this->rasterBandRead = std::vector<bool>(this->getBandCount(), true);

	//initialize display raster band properties
	this->displayRasterBandPointers = std::vector<void *>(this->getBandCount(), nullptr);
	this->displayRasterBandRead = std::vector<bool>(this->getBandCount(), false);
}

/******************************************************************************
			      ~GDALRasterWrapper()
******************************************************************************/
GDALRasterWrapper::~GDALRasterWrapper() {
	if (this->rasterAllocated) {
		CPLFree(p_raster);
	}

	if (this->displayRasterAllocated) {
		CPLFree(p_displayRaster);
	}
}

/******************************************************************************
				  getDataset()
******************************************************************************/
GDALDataset *GDALRasterWrapper::getDataset() {
	return this->p_dataset.get();
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
				   getBandNames()
******************************************************************************/
std::vector<std::string> GDALRasterWrapper::getBands() {
	std::vector<std::string> retval;

	for( auto&& p_band : this->p_dataset->GetBands() ){
		retval.push_back(p_band->GetDescription());
	}

	return retval;
}

/******************************************************************************
			       getGeotransform()
******************************************************************************/
double *GDALRasterWrapper::getGeotransform() {
	return this->geotransform;
}

/******************************************************************************
				getMinPixelVal()
******************************************************************************/
double GDALRasterWrapper::getMinPixelVal(int band) {
	return this->p_dataset->GetRasterBand(band + 1)->GetMinimum();
}

/******************************************************************************
				getMaxPixelVal()
******************************************************************************/
double GDALRasterWrapper::getMaxPixelVal(int band) {
	return this->p_dataset->GetRasterBand(band + 1)->GetMaximum();
}

/******************************************************************************
				allocateRaster()
******************************************************************************/
void GDALRasterWrapper::allocateRaster(bool display) {
	//This is not an ideal usage of CPLMalloc because it may be very large
	//In the future, this will likely be swapped out.
	if (display) {
		size_t bandSize = this->getRasterTypeSize() * this->displayRasterWidth * this->displayRasterHeight;
		this->p_displayRaster = CPLRealloc(this->p_displayRaster, bandSize * this->getBandCount());
		this->displayRasterAllocated = true;
		for (int i = 0; i < this->getBandCount(); i++) {
			this->displayRasterBandPointers[i] = (void *)((size_t)this->p_displayRaster + ((size_t)i * bandSize));
			this->displayRasterBandRead[i] = false;
		}
	}
	else {
		size_t bandSize = this->getRasterTypeSize() * this->getWidth() * this->getHeight();
		this->p_raster = CPLMalloc(bandSize * this->getBandCount());
		this->rasterAllocated = true;
		for (int i = 0; i < this->getBandCount(); i++) {
			this->rasterBandPointers[i] = (void *)((size_t)this->p_raster + ((size_t)i * bandSize));
		}
	}
}

/******************************************************************************
			    readRasterBand()
******************************************************************************/
void GDALRasterWrapper::readRasterBand(void *p_band, int width, int height, int band) {	
	/**
	 * see https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec
	 * and https://gdal.org/en/stable/tutorials/raster_api_tut.html#reading-raster-data
	 */

	//perform raster read on current band
	CPLErr err = this->p_dataset->GetRasterBand(band + 1)->RasterIO(
		GF_Read, 			//GDALRWFlag eRWFlag
		0, 				//int nXOff
		0,				//int nYOff
		this->getWidth(),		//int nXSize
		this->getHeight(),		//int nYSize
		p_band,				//void *pData
		width,				//int nBufXSize
		height,				//int nBufYSize
		this->getRasterType(),		//GDALDataType eBufType 
		0,				//int nPixelSpace
		0				//int nLineSpace
	);
	if (err) {
		throw std::runtime_error("error reading raster band from dataset.");
	}
}

/******************************************************************************
				  getBuffer()
******************************************************************************/
template <typename T> 
py::buffer GDALRasterWrapper::getBuffer(size_t size, void *p_buffer, int width, int height) {
	//see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	return py::memoryview::from_buffer(
		(T*)p_buffer, 					//buffer
		{this->getBandCount(), height, width}, 		//shape
		{size * height * width, size * width, size} 	//stride
	);
}	

/******************************************************************************
			      getRasterAsMemView()
******************************************************************************/
py::buffer GDALRasterWrapper::getRasterAsMemView(int width, int height) {
	bool display = (width != this->getWidth() || height != this->getHeight());
	void *p_buffer;

	//allocate raster if required
	if (!display && !this->rasterAllocated) {
		this->allocateRaster(display); //allocate whole raster
	}

	//reallocate display raster if required
	if (display && 
		(!this->displayRasterAllocated || 
		 width != this->displayRasterWidth || 
		 height != this->displayRasterHeight)) 
	{
		this->displayRasterWidth = width;
		this->displayRasterHeight = height;
		this->allocateRaster(display); //reallocate (potentially downsampled) display raster
	}
	
	//read raster bands if required
	if (!display) {
		for (int i = 0; i < this->getBandCount(); i++) {
			if (!this->rasterBandRead[i]) {
				this->readRasterBand(this->rasterBandPointers[i], width, height, i);
				this->rasterBandRead[i] = true;
			}
		}
		p_buffer = this->p_raster;
	}
	else {
		for (int i = 0; i < this->getBandCount(); i++) {
			if (!this->displayRasterBandRead[i]) {
				this->readRasterBand(this->displayRasterBandPointers[i], width, height, i);
				this->displayRasterBandRead[i] = true;
			}
		}
		p_buffer = this->p_displayRaster;
	}
	
	switch(this->getRasterType()) {
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

/******************************************************************************
				getRasterBand()				     
******************************************************************************/
void *GDALRasterWrapper::getRasterBand(int band) {
	if (!this->rasterAllocated) {
		this->allocateRaster(false); //display = false
	}

	if (!this->rasterBandRead[band]) {
		this->readRasterBand(
			this->rasterBandPointers[band], 
			this->getWidth(),
			this->getHeight(),
			band
		);
		this->rasterBandRead[band] = true;
	}

	return this->rasterBandPointers[band];
}

/******************************************************************************
				getRasterType()				     
******************************************************************************/
GDALDataType GDALRasterWrapper::getRasterType() { 
	return this->rasterType.GetNumericDataType(); 
}

/******************************************************************************
			      getRasterTypeSize()				     
******************************************************************************/
size_t GDALRasterWrapper::getRasterTypeSize() {
	return this->rasterType.GetSize();
}

/******************************************************************************
			      getMinIndexIntType()		     
******************************************************************************/
std::string GDALRasterWrapper::getMinIndexIntType(bool singleBand) {
	size_t maxIndex = this->getWidth() * this->getHeight();
	maxIndex = singleBand ? maxIndex : maxIndex * this->getBandCount();
	if (std::numeric_limits<unsigned short>::max() >= maxIndex) {
		return "unsigned_short";
	}
	else if (std::numeric_limits<unsigned>::max() >= maxIndex) {
		return "unsigned";
	}
	else if (std::numeric_limits<unsigned long>::max() >= maxIndex) {
		return "unsigned_long";
	}
	else {
		return "unsigned_long_long";
	}
}

/******************************************************************************
				    write()				     
******************************************************************************/
void GDALRasterWrapper::write(std::string filename) {
	std::filesystem::path filepath = filename;
	std::string extension = filepath.extension().string();

	//TODO move this to write.cpp ? a write_raster() function ??
	if (extension != ".tif") {
		throw std::runtime_error("write only supports .tif files right now");
	}

	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("GTiff");

	GDALClose(p_driver->CreateCopy(filename.c_str(), this->p_dataset.get(), (int)false, nullptr, nullptr, nullptr));
}
