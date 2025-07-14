/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of systematic sampling
 * Author: Joseph Meyer
 * Date: July, 2025

 *
 ******************************************************************************/

#include <iostream> //TODO remove
#include <random>

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
	
	//TODO generate random corner (origin) for the grid
	double xDiff = xMax - xMin;
	double yDiff = yMax - yMin;
	double rngMax = yDiff * xDiff;
	std::mt19937::result_type seed = time(nullptr);
	auto rng = std::bind(
		std::uniform_real_distribution<double>(0, rngMax),
		std::mt19937(seed)
	);

	double yCoord = rng() / xDiff; //divide by xDiff because the result will then be between 0 and yDiff
	double xCoord = rng() / yDiff; //divide by yDiff becausethe result will then be between 0 and xDiff

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

	//query to create grid
	OGRLayer *p_gridTest = p_raster->getDataset()->ExecuteSQL(queryString.c_str(), nullptr, "SQLITE");

	//wktPoints represents the samples as well known text, and is returned to the (Python) caller
	std::vector<std::string> wktPoints;	

	//coordinate representations the samples as vectors only if PLOT is true, returned to the (Python) caller
	std::vector<double> xCoords, yCoords;

	//grid represents the grid used to create the sampels only if PLOT is true, and is returned to the (Python) caller
	std::vector<std::vector<std::vector<double>>> grid;

	//iterate through the polygons in the grid and populate the ruturn data depending on user inputs
	for (const auto& p_feature : *p_gridTest) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		for (const auto& p_polygon : *p_geometry->toMultiPolygon()) {
			//generate sample depending on 'location' parameter
			OGRPoint point;
			if (location == "centers") {
				p_polygon->Centroid(&point);
			}
			else if (location == "corners") {
				(*p_polygon->begin())->StartPoint(&point);
			}
			else { //location == "random"
				OGREnvelope envelope;
				p_polygon->getEnvelope(&envelope);
				double xMinEnv = envelope.MinX;
				double xMaxEnv = envelope.MaxX;
				double xDiffEnv = xMaxEnv - xMinEnv;
				double yMinEnv = envelope.MinY;
				double yMaxEnv = envelope.MaxY;
				double yDiffEnv = yMaxEnv - yMinEnv;

				point.setX(xMinEnv + rng() / (rngMax / xDiffEnv));
				point.setY(yMinEnv + rng() / (rngMax / yDiffEnv));
				while (!p_polygon->Contains(&point)) {
					point.setX(xMinEnv + rng() / (rngMax / xDiffEnv));
					point.setY(yMinEnv + rng() / (rngMax / yDiffEnv));
				}
			}
			double x = point.getX();
			double y = point.getY();

			//if point is not within raster extent, don't add it
			//TODO don't necessarily just continue here, because we may have to add grid to print
			if (x < xMin || x > xMax || y < yMin || y > yMax) {
				continue;
			}

			wktPoints.push_back(point.exportToWkt());

			//if we're plotting, add point to coordinate vectors and plottable vectors to 'grid'
			if (plot) {
				//for point plotting
				xCoords.push_back(x);
				yCoords.push_back(y);

				//for grid plotting
				grid.push_back({}); //add new polygon to grid
				grid.back().push_back({}); //add new x vector to polygon grid
				grid.back().push_back({}); //add new y vector to polygon grid
				for (const auto& p_linearRing : *p_polygon) {
					for (const auto& p_point : *p_linearRing) {
						grid.back()[0].push_back(p_point.getX());
						grid.back()[1].push_back(p_point.getY());
					}
				}
			}
		}
	}

	return {wktPoints, {xCoords, yCoords}, grid};
}

PYBIND11_MODULE(systematic, m) {
	m.def("systematic_cpp", &systematic);
}
