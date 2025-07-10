/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of systematic sampling
 * Author: Joseph Meyer
 * Date: July, 2025

 *
 ******************************************************************************/

#include <iostream> //TODO remove

#include "access.h"
#include "raster.h"
#include "vector.h"
#include "write.h"

/*
 *
 */
std::tuple<
	std::vector<std::string>, //wkt samples
	std::vector<std::vector<double>>, //array of samples to plot
	std::vector<std::vector<std::vector<double>>> //array of grid to plot
>
systematic(
	GDALRasterWrapper *p_raster,
	double cellsize,
	std::string shape,
	std::string location,
	bool plot,
	std::string filename)
{
	//Step 1: formulate sql query
	double size = 200; //use cellsize for real query
	//calculate extent polygon
	double xMin, xMax, yMin, yMax;
	xMin = p_raster->getXMin();
	xMax = p_raster->getXMax();
	yMin = p_raster->getYMin();
	yMax = p_raster->getYMax();
	std::string extentPolygon = "'POLYGON (( " 
		+ std::to_string(xMin) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMin) + " ))'";

	std::string queryString = "SELECT ST_SquareGrid(ST_GeomFromTest("
		+ extentPolygon
		+ "), " + std::to_string(size) + ")";

	std::cout << "query string: " << std::endl;
	std::cout << queryString << std::endl;

	OGRLayer *p_gridTest = p_raster->getDataset()->ExecuteSQL(queryString.c_str(), nullptr, "SQLITE");

	std::cout << "HERE!" << std::endl;
	std::cout << p_gridTest << std::endl;
	
	// SELECT ST_HexagonalGrid(ST_GeomFromText(), 1)	

	//step 2: use query to generate points
	return {{""},{{0.0}},{{{0.0}}}};
}

PYBIND11_MODULE(systematic, m) {
	m.def("systematic_cpp", &systematic);
}
