#include "raster.h"
#include <iostream>


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
}

std::string GDALRasterWrapper::getDriver() {
	return std::string(this->p_dataset->GetDriverName());
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

double GDALRasterWrapper::getPixelHeight() {
	return std::abs(this->geotransform[1]);
}

double GDALRasterWrapper::getPixelWidth() {
	return std::abs(this->geotransform[5]);
}

std::vector<std::string> GDALRasterWrapper::getBands() {
	std::vector<std::string> retval;

	for( auto&& p_band : this->p_dataset->GetBands() ){
		retval.push_back(p_band->GetDescription());
	}

	return retval;
}
