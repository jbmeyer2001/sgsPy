#pragma once

#include <string>
#include <vector>
#include <unordered_map>

#include <Python.h>
#include <gdal.h>

class Raster {
	GDALDataset* dataset;
	double geotransform[6];
	std::string driver;
	std::string crs;
	int width;
	int height;
	int layers;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double pixelHeight;
	double pixelWidth;
	std::vector<std::string> bands;
	std::unordered_map<std::string, int> bandNameMap;

	//virtual memory will be used (for when there are large files)
	//https://github.com/OSGeo/gdal/blob/5ebd51a11c939b1b2a62002e25127beddd2f1da0/gcore/gdal.h#L2339
	//https://gdal.org/en/latest/doxygen/cpl__virtualmem_8h.html
	
	Raster(std:::string
};
