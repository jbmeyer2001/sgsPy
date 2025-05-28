#include <iostream>
#include <stdexcept>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

py::module_ json = py::module_::import("json");

/**
 * 
 */
class SpatialRaster {
	private:
	std::unique_ptr<GDALDataset> p_dataset;
	double geotransform[6];
	std::string driver;
	std::string crs;
	int width;
	int height;
	int layers;
	double xMin;
	double xMax;
	double yMin;
	double yMax;
	double pixelHeight;
	double pixelWidth;
	std::vector<std::string> bands;
	std::unordered_map<std::string, int> bandNameMap;

	public:

	/**
	 *
	 */
	SpatialRaster(std::string filename) {
		//dataset
		this->p_dataset = std::unique_ptr<GDALDataset>(GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly)));
		
		//geotransfrom
		CPLErr cplerr = this->p_dataset->GetGeoTransform(this->geotransform);
		if (cplerr) {
			throw std::runtime_error("error getting geotransform from dataset.");
		}

		//metadata
		this->driver = std::string(p_dataset->GetDriverName());
		char *p_crs = this->crs.data();
		OGRErr ogrerr = OGRSpatialReference(this->p_dataset->GetProjectionRef()).exportToPROJJSON(&p_crs, nullptr);
		if (ogrerr) {
			throw std::runtime_error("error getting coordinate reference system from dataset.");
		}
		this->crs = std::string(p_crs);
			
		//raster(s) dimensions
		this->width = p_dataset->GetRasterXSize();
		this->height = p_dataset->GetRasterYSize();
		this->layers = p_dataset->GetRasterCount();

		//raster extent
		double xbounds[2] = {
			this->geotransform[0], 
			this->geotransform[0] + this->geotransform[1] * this->width + this->geotransform[2] * this->height
		};
		double ybounds[2] = {
			this->geotransform[3],
			this->geotransform[3] + this->geotransform[4] * this->width + this->geotransform[5] * this->height
		};
		this->xMin = std::min(xbounds[0], xbounds[1]);
		this->xMax = std::max(xbounds[0], xbounds[1]);
		this->yMin = std::min(ybounds[0], ybounds[1]);
		this->yMax = std::max(ybounds[0], ybounds[1]);
		this->pixelHeight = std::abs(this->geotransform[1]);
		this->pixelWidth = std::abs(this->geotransform[5]);

		//bands
		for (int i; i < this->layers; i++) {
			std::string bandName = std::string(this->p_dataset->GetRasterBand(i)->GetDescription());
			this->bands.push_back(bandName);
			this->bandNameMap.emplace(bandName, i);
		}
	}

	/**
	 *
	 */
	std::string getDriver() { 
		return this->driver; 
	}

	/**
	 *
	 */
	py::dict getCRS() { return json.attr("loads")(this->crs); }
	
	/**
	 *
	 */
	int getHeight() { return this->height; }
	
	/**
	 *
	 */
	int getWidth() { return this->width; }
	
	/**
	 *
	 */
	int getLayers() { return this->layers; }
	
	/**
	 *
	 */
	double getXMin() { return this->xMin; }
	
	/**
	 *
	 */
	double getXMax() { return this->xMax; }
	
	/**
	 *
	 */
	double getYMin() { return this->yMin; }
	
	/**
	 *
	 */
	double getYMax() { return this->yMax; }
	
	/**
	 *
	 */
	double getPixelHeight() { return this->pixelHeight; }
	
	/**
	 *
	 */
	double getPixelWidth() { return this->pixelWidth; }
	
	/**
	 *
	 */
	std::vector<std::string> getBands() { return this->bands; }

	/*
	template <typename T>
	//TODO find some way to make this accessable by the python [] operator
	std::vector<T> getRasterAsNumpy() {
		//TODO add
		return std::vector<int>;
	}

	std::vector<T> getVirtualMemoryRaster() {
		//TODO add
		return std::vector<int>;
	}
	*/
};


/**
 *
 */
PYBIND11_MODULE(SpatialRaster, m) {
	py::class_<SpatialRaster>(m, "SpatialRaster")
		.def(py::init<std::string>())
		.def("driver", &SpatialRaster::getDriver)
		.def("crs", &SpatialRaster::getCRS)
		.def("height", &SpatialRaster::getHeight)
		.def("width", &SpatialRaster::getWidth)
		.def("layers", &SpatialRaster::getLayers)
		.def("xmin", &SpatialRaster::getXMin)
		.def("xmax", &SpatialRaster::getXMax)
		.def("ymin", &SpatialRaster::getYMin)
		.def("ymax", &SpatialRaster::getYMax)
		.def("pixel_height", &SpatialRaster::getPixelHeight)
		.def("pixel_width", &SpatialRaster::getPixelWidth)
		.def("bands", &SpatialRaster::getBands);
		//.def( something for getting the raster as a numpy array
}
