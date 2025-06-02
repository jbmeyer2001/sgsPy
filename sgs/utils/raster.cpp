#include "raster.h"

GDALRasterWrapper::GDALRasterWrapper(std::string filename) {
	//must register drivers first
	GDALAllRegister();

	//dataset
	p_dataset = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly)));
	
	//geotransform
	CPLErr cplerr = this->p_dataset->GetGeoTransform(this->geotransform);
	if (cplerr) {
		throw std::runtime_error("error getting geotransform from dataset");
	}	

	//gdal data type
	this->type = this->p_dataset->GetRasterBand(1)->GetRasterDataType();	
}

GDALRasterWrapper::~GDALRasterWrapper() {
	if (this->p_raster != nullptr) {
		CPLFree(this->p_raster);
	}
}

std::string GDALRasterWrapper::getDriver() {
	return std::string(this->p_dataset->GetDriverName()) 
		+ "/" 
		+ std::string(this->p_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
}

std::string GDALRasterWrapper::getCRS() {
	char *p_crs;

	OGRErr ogrerr = OGRSpatialReference(this->p_dataset->GetProjectionRef()).exportToPROJJSON(&p_crs, nullptr);
	if (ogrerr) {
		throw std::runtime_error("error getting coordinate reference system from dataset.");
	}
		
	return std::string(p_crs);
}

int GDALRasterWrapper::getWidth() {
	return this->p_dataset->GetRasterXSize();
}

int GDALRasterWrapper::getHeight() {
	return this->p_dataset->GetRasterYSize();
}

int GDALRasterWrapper::getLayers() {
	return this->p_dataset->GetRasterCount();
}

double GDALRasterWrapper::getXMax() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::max(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
	);
}

double GDALRasterWrapper::getXMin() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::min(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
	);
}

double GDALRasterWrapper::getYMax() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::max(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
	);
}

double GDALRasterWrapper::getYMin() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::min(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
	);
}

double GDALRasterWrapper::getPixelWidth() {
	return std::abs(this->geotransform[5]);
}

double GDALRasterWrapper::getPixelHeight() {
	return std::abs(this->geotransform[1]);
}

std::vector<std::string> GDALRasterWrapper::getBands() {
	std::vector<std::string> retval;

	for( auto&& p_band : this->p_dataset->GetBands() ){
		retval.push_back(p_band->GetDescription());
	}

	return retval;
}

void GDALRasterWrapper::allocateRasterHelper(size_t size) {
	//allocate the whole raster
	this->p_raster = CPLMalloc(size * this->getLayers() * this->getHeight() * this->getWidth());

	//calculate the number of bytes that a single raster layer/band takes
	size_t layerSize = size * this->getHeight() * this->getWidth();

	//iterate and read the raster bands
	//NOTE: GDALRasterBand::RasterIO takes care of windowing for us so we *should* be able
	//to allocate a raster larger than the size of RAM... (should)
	void *bufferLocation = this->p_raster;
	for( auto&& p_band : this->p_dataset->GetBands() ) {
		p_band->RasterIO(
			GF_Read, 			//GDALRWFlag eRWFlag
			0, 				//int nXOff
			0,				//int nYOff
			this->getWidth(),		//int nXSize
			this->getHeight(),		//int nYSize
			bufferLocation,			//void *pData
			this->getWidth(),		//int nBufXSize
			this->getHeight(),		//int nBufYSize
			this->type,			//GDALDataType eBufType 
			0,				//int nPixelSpace
			0				//int nLineSpace
		);
		bufferLocation = (void *)((size_t)bufferLocation + layerSize);
	}

}

void GDALRasterWrapper::allocateRaster() {
	if (this->rasterAllocated) {
		throw std::runtime_error("cannot allocate an already allocated raster");
	}
	
	switch(this->type) {
		case GDT_Int8:
			this->allocateRasterHelper(1);
			break;
		case GDT_UInt16:
		case GDT_Int16:
			this->allocateRasterHelper(2);
			break;
		case GDT_UInt32:
		case GDT_Int32:
		case GDT_Float32:
			this->allocateRasterHelper(4);
			break;
		case GDT_UInt64:
		case GDT_Int64:
		case GDT_Float64:
			this->allocateRasterHelper(8);
			break;
		default:
			throw std::runtime_error("raster pixel data type not acceptable.");
	}

	this->rasterAllocated = true;
}

template <typename T> py::buffer GDALRasterWrapper::getBuffer(size_t size) {
	//see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	return py::memoryview::from_buffer(
		(T*)this->p_raster, //buffer
		{this->getLayers(), this->getHeight(), this->getWidth()}, //shape
		{size * this->getHeight() * this->getWidth(), size * this->getWidth(), size} //stride
	);
}	

py::buffer GDALRasterWrapper::getRasterAsMemView() {
	if (!this->rasterAllocated) {
		this->allocateRaster();
	}
	
	switch(this->type) {
		case GDT_Int8:
			return getBuffer<int8_t>(1);
		case GDT_UInt16:
			return getBuffer<uint16_t>(2);
		case GDT_Int16:
			return getBuffer<int16_t>(2);
		case GDT_UInt32:
			return getBuffer<uint32_t>(4);
		case GDT_Int32:
			return getBuffer<int32_t>(4);
		case GDT_Float32:
			return getBuffer<float>(4);
		case GDT_UInt64:
			return getBuffer<uint64_t>(8);
		case GDT_Int64:
			return getBuffer<int64_t>(8);
		case GDT_Float64:
			return getBuffer<double>(8);
		default:
			throw std::runtime_error("raster pixel data type not acceptable.");
	}
}

void *GDALRasterWrapper::getRaster() {
	if (!this->rasterAllocated) {
		this->allocateRaster();
	}

	return this->p_raster;
}
