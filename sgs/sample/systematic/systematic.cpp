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
	double cellSize,
	std::string shape,
	std::string location,
	bool plot,
	std::string filename)
{

	
	//Step 1: formulate sql query
	//calculate extent polygon
	double xMin, xMax, yMin, yMax;
	xMin = p_raster->getXMin();
	xMax = p_raster->getXMax();
	yMin = p_raster->getYMin();
	yMax = p_raster->getYMax();
	
	//determine grid creation function
	std::string gridFunction;
	if (shape == "square") {
		gridFunction = "ST_SquareGrid";
	}
	else if (shape == "hexagon") {
		gridFunction = "ST_HexagonalGrid";
	}
	else { //shape == "triangle"
		gridFunction = "ST_TriangularGrid";
	}

	//create sql query using extent polygon and grid function
	std::string extentPolygon = "'POLYGON (( " 
		+ std::to_string(xMin) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMin) + " ))'";

	std::string queryString = "SELECT " + gridFunction + "(ST_GeomFromText("
		+ extentPolygon
		+ "), " + std::to_string(cellSize) + ")";

	std::cout << "query string: " << std::endl;
	std::cout << queryString << std::endl << std::endl;
	OGRLayer *p_gridTest = p_raster->getDataset()->ExecuteSQL(queryString.c_str(), nullptr, "SQLITE");
	std::cout << p_gridTest << std::endl;
	std::cout << "HERE!" << std::ends;

	std::vector<std::vector<std::vector<double>>> grid;
	for (const auto& p_feature : *p_gridTest) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		std::cout << p_geometry->getGeometryName() << std::endl;
		for (const auto& p_polygon : *p_geometry->toMultiPolygon()) {
			for (const auto& p_linearRing : *p_polygon) {
				std::vector<double> xCoords;
				std::vector<double> yCoords;
				for (const auto& p_point : *p_linearRing) {
					xCoords.push_back(p_point.getX());
					yCoords.push_back(p_point.getY());
				}
				grid.push_back({xCoords, yCoords});
			}
		}
	}

	//step 2: use query to generate points
	return {{""},{{0.0}},grid};
}

PYBIND11_MODULE(systematic, m) {
	m.def("systematic_cpp", &systematic);
}
