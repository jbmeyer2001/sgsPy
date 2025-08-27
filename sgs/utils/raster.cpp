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
	GDALDataset *p_dataset = GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly));
	if (!p_dataset) {
		throw std::runtime_error("dataset pointer is null after initialization, dataset unable to be initialized.");
	}

	this->createFromDataset(p_dataset);
}	

/******************************************************************************
			      GDALRasterWrapper()
******************************************************************************/
GDALRasterWrapper::GDALRasterWrapper(GDALDataset *p_dataset, std::vector<void *>bands) {
	this->createFromDataset(p_dataset);
	this->rasterBandPointers = bands;
	this->rasterBandRead = std::vector<bool>(bands.size(), true);
}

/******************************************************************************
			      GDALRasterWrapper()
******************************************************************************/
GDALRasterWrapper::GDALRasterWrapper(GDALDataset *p_dataset) {
	this->createFromDataset(p_dataset);
}


/******************************************************************************
			      createFromDataset()
******************************************************************************/
void GDALRasterWrapper::createFromDataset(GDALDataset *p_dataset) {
	this->p_dataset = GDALDatasetUniquePtr(p_dataset);

	//geotransform
	CPLErr cplerr = this->p_dataset->GetGeoTransform(this->geotransform);
	if (cplerr) {
		throw std::runtime_error("error getting geotransform from dataset.");
	}

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
	std::vector<GDALDataType> rasterBandTypes,
	int width,
	int height,
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
		rasterBandTypes[0],
		nullptr
	));
	if (!this->p_dataset) {
		throw std::runtime_error("dataset pointer is null after initialization, dataset unable to be initialized.");
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
		err = this->p_dataset->AddBand(rasterBandTypes[i], papszOptions);
		if (err) {
			throw std::runtime_error("error adding band.");
		}
		GDALRasterBand *p_band = this->p_dataset->GetRasterBand(i + 1);
		p_band->SetDescription(rasterBandNames[i].c_str());
		p_band->SetNoDataValue(-1); //all strat rasters processed by this package will have nodata -1 
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
	
	//initialize raster band properties
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
	for (int i = 0; i < this->getBandCount(); i++) {
		if (this->rasterBandRead[i]) {
			CPLFree(this->rasterBandPointers[i]);
		}

		if (this->displayRasterBandRead[i]) {
			CPLFree(this->displayRasterBandPointers[i]);
		}
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
			       getBandNoDataValue()
******************************************************************************/
double GDALRasterWrapper::getBandNoDataValue(int band) {
	GDALRasterBand *p_band = this->getRasterBand(band);
	return p_band->GetNoDataValue();
}

/******************************************************************************
			    readRasterBand()
******************************************************************************/
void GDALRasterWrapper::readRasterBand(int width, int height, int band) {	
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

/******************************************************************************
				  getBuffer()
******************************************************************************/
template <typename T> 
py::buffer GDALRasterWrapper::getBuffer(size_t size, void *p_buffer, int width, int height) {
	//see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	return py::memoryview::from_buffer(
		(T*)p_buffer, 		//buffer
		{height, width}, 	//shape
		{size * width, size} 	//stride
	);
}	

/******************************************************************************
			    getRasterBandAsMemView()
******************************************************************************/
py::buffer GDALRasterWrapper::getRasterBandAsMemView(int width, int height, int band) {
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

/******************************************************************************
				getRasterBand()				     
******************************************************************************/
GDALRasterBand *GDALRasterWrapper::getRasterBand(int band) {
	return this->p_dataset->GetRasterBand(band + 1);	
}

/******************************************************************************
			     getRasterBandBuffer()				     
******************************************************************************/
void *GDALRasterWrapper::getRasterBandBuffer(int band) {
	if (!this->rasterBandRead[band]) {
		this->readRasterBand(
			this->getWidth(), 
			this->getHeight(),
			band
		);
	}

	return this->rasterBandPointers[band];
}

/******************************************************************************
			      getRasterBandType()				     
******************************************************************************/
GDALDataType GDALRasterWrapper::getRasterBandType(int band) {
	GDALRasterBand *p_band = this->p_dataset->GetRasterBand(band + 1);
	return p_band->GetRasterDataType();	
}

/******************************************************************************
			    getRasterBandTypeSize()				     
******************************************************************************/
size_t GDALRasterWrapper::getRasterBandTypeSize(int band) {
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

/******************************************************************************
				    write()				     
******************************************************************************/
void GDALRasterWrapper::write(std::string filename) {
	std::filesystem::path filepath = filename;
	std::string extension = filepath.extension().string();

	if (extension != ".tif") {
		throw std::runtime_error("write only supports .tif files right now");
	}

	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("GTiff");

	GDALClose(p_driver->CreateCopy(filename.c_str(), this->p_dataset.get(), (int)false, nullptr, nullptr, nullptr));
}
