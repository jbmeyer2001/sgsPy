#include "raster.h"

py::module_ json = py::module_::import("json");

SpatialRaster::SpatialRaster(std::string filename) {
	//dataset
	this->p_dataset = std::unique_ptr<GDALDataset>(GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly)));
		
	//geotransfrom
	CPLErr cplerr = this->p_dataset->GetGeoTransform(this->geotransform);
	if (cplerr) {
		throw std::runtime_error("error getting geotransform from dataset.");
	}
}

std::string SpatialRaster::getDriver() {
	return std::string(this->p_dataset->GetDriverName());
}

py::dict SpatialRaster::getCRS() {
	char *p_crs;

	OGRErr ogrerr = OGRSpatialReference(this->p_dataset->GetProjectionRef()).exportToPROJJSON(&p_crs, nullptr);
	if (ogrerr) {
		throw std::runtime_error("error getting coordinate reference system from dataset.");
	}
		
	return json.attr("loads")(std::string(p_crs));
}

int SpatialRaster::getWidth() {
	return this->p_dataset->GetRasterXSize();
}

int SpatialRaster::getHeight() {
	return this->p_dataset->GetRasterYSize();
}

int SpatialRaster::getLayers() {
	return this->p_dataset->GetRasterCount();
}

double SpatialRaster::getXMax() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::max(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
	);
}

double SpatialRaster::getXMin() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::min(
		this->geotransform[0], 
		this->geotransform[0] + this->geotransform[1] * width + this->geotransform[2] * height
	);
}

double SpatialRaster::getYMax() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::max(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
	);
}

double SpatialRaster::getYMin() {
	int width = this->p_dataset->GetRasterXSize();
	int height = this->p_dataset->GetRasterYSize();
	return std::min(
		this->geotransform[3],
		this->geotransform[3] + this->geotransform[4] * width + this->geotransform[5] * height
	);
}

double SpatialRaster::getPixelHeight() {
	return std::abs(this->geotransform[1]);
}

double SpatialRaster::getPixelWidth() {
	return std::abs(this->geotransform[5]);
}

std::vector<std::string> SpatialRaster::getBands() {
	std::vector<std::string> retval;

	for( auto&& p_band : this->p_dataset->GetBands() ){
		retval.push_back(p_band->GetDescription());
	}

	return retval;
}
