#include "raster.h"

py::module_ json = py::module_::import("json");

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
	if (this->p_CPLVirtualMemRaster != nullptr) {
		CPLVirtualMemFree(this->p_CPLVirtualMemRaster);
	}
}

std::string GDALRasterWrapper::getDriver() {
	return std::string(this->p_dataset->GetDriverName()) 
		+ "/" 
		+ std::string(this->p_dataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
}

py::dict GDALRasterWrapper::getCRS() {
	char *p_crs;

	OGRErr ogrerr = OGRSpatialReference(this->p_dataset->GetProjectionRef()).exportToPROJJSON(&p_crs, nullptr);
	if (ogrerr) {
		throw std::runtime_error("error getting coordinate reference system from dataset.");
	}
		
	return json.attr("loads")(std::string(p_crs));
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

void GDALRasterWrapper::allocateRaster() {
	if (this->p_CPLVirtualMemRaster != nullptr) {
		throw std::runtime_error("raster already allocated, cannot call allocateRaster() again.");
	}

	//see https://github.com/OSGeo/gdal/blob/9f4bc2d28f853d9e39a59656cb8f3318b51f9be2/gcore/gdalvirtualmem.cpp#L764
	//TODO dynamically get physical memory size to pass to nCacheSize parameter
	this->p_CPLVirtualMemRaster = GDALDatasetGetVirtualMem(this->p_dataset.get(), 					//GDALDatasetH hDS
		GF_Read,						//GDALRWFlag eRWFlag
		0,									//int nXOff
		0,									//int nYOff
		this->p_dataset->GetRasterXSize(),	//int nXSize
		this->p_dataset->GetRasterYSize(),	//int nYSize
		this->p_dataset->GetRasterXSize(),	//int nBufXSize
		this->p_dataset->GetRasterYSize(),	//int nBufYSize
		this->type,							//GDALDataType eBufType
		this->p_dataset->GetRasterCount(),	//int nBandCount
		NULL,								//int *panBandMap
		0,									//int nPixelSpace	
		0,									//GIntBig nLineSpace
		0,									//GIntBig nBandSpace
		10485760, 							//size_t nCacheSize
		0,									//size_t nPageSizeHint
		(int)false,							//int bSingleThreadUsage
		NULL);								//CSLConstList papszOptions
}

void *GDALRasterWrapper::getRasterPointer() {
	if (this->p_CPLVirtualMemRaster == nullptr) {
		this->allocateRaster();
	}

	return CPLVirtualMemGetAddr(this->p_CPLVirtualMemRaster);
}

template <typename T> py::buffer GDALRasterWrapper::getBuffer(size_t size) {
	//see https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html#memory-view
	return py::memoryview::from_buffer(
		(T*)this->getRasterPointer(), //buffer
		{this->getLayers(), this->getHeight(), this->getWidth()}, //shape
		{size * this->getHeight() * this->getWidth(), size * this->getWidth(), size} //stride
	);
}	

py::buffer GDALRasterWrapper::getRasterAsMemView() {
	size_t size;
	std::string formatDescriptor;

	//data size (bytes) and corresponding python type depends on the gdal data type
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
		return this->getRasterPointer();
}
