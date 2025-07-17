/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of systematic sampling
 * Author: Joseph Meyer
 * Date: July, 2025

 *
 ******************************************************************************/

#include <iostream>
#include <random>

#include "raster.h"
#include "vector.h"
#include "write.h"

/*
 * This function conducts Systematic sampling on an input raster image.
 *
 * First, the extent polygon of the raster is determined, then a random
 * point is found within the polygon to act as the origin. An SQL query
 * is conducted using one of the ST_SquareGrid, or ST_HexagonalGrid
 * spatialite functions to create a grid of polygons
 * of the user-specified shape. The grid is then rotated by a randomly
 * generated rotation angle.
 *
 * Next, the resulting grid polygons are iterated through, and a sample
 * point is determined for each polygon depending on the user-defined
 * location parameter (centers, corners, or random). Plot-required data
 * is saved if plot is true (to later be utilized by the Python side
 * of the application with matplotlib). Additionally, a vector of 
 * points may be generated and used to write if the user specified a
 * valid filename.
 *
 * @param GDALRasterWrapper *p_raster raster to be systematically sampled
 * @param double cellSize the size of the grid cell shapes
 * @param std::string shape the shape of the grid cells
 * @param std::string location the location within a cell to sample
 * @param bool plot whether to save and return plot-required data
 * @param std::string filename to write to or "" if not to write
 * @returns std::tuple<
 * 		std::vector<std::string>,
 * 		std::vector<std::vector<double>>,
 * 		std::vector<std::vector<std::vector<double>>>
 * 	>
 * 	the vector of strings is a wkt strings of the sample points,
 * 	the 2d vector of doubles contains sample points to plot
 * 	the 3d vector of doubles contains grid cells to plot
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
	
	//determine raster extent
	double xMin, xMax, yMin, yMax;
	xMin = p_raster->getXMin();
	xMax = p_raster->getXMax();
	yMin = p_raster->getYMin();
	yMax = p_raster->getYMax();
	
	//generate random number generator
	double xDiff = xMax - xMin;
	double yDiff = yMax - yMin;
	double rngMax = yDiff * xDiff;
	std::mt19937::result_type seed = time(nullptr);
	auto rng = std::bind(
		std::uniform_real_distribution<double>(0, rngMax),
		std::mt19937(seed)
	);

	//determine random origin location within extent polygon
	double yCoord = rng() / xDiff; //divide by xDiff because the result will then be between 0 and yDiff
	double xCoord = rng() / yDiff; //divide by yDiff becausethe result will then be between 0 and xDiff

	//determine random rotation angle within extent polygon
	double rotation = rng() / (rngMax / 180);

	//determine grid creation function
	std::string gridFunction;
	if (shape == "square") {
		gridFunction = "ST_SquareGrid";
	}
	else { //shape == "hexagon" 
		gridFunction = "ST_HexagonalGrid";
	}

	//create sql query using extent polygon and grid function
	std::string extentPolygon = "'POLYGON (( " 
		+ std::to_string(xMin) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMin) + " ))'";

	std::string queryString = "SELECT RotateCoords(" 
		+ gridFunction 
		+ "(RotateCoords("
		+ "ST_GeomFromText("
		+ extentPolygon
		+ "), " + std::to_string(-rotation) + "), "
		+ std::to_string(cellSize) + "), "
		+ std::to_string(rotation) + ")";

	std::cout << "query string:" << std::endl;
	std::cout << queryString << std::endl;	

	//query to create grid
	OGRLayer *p_gridTest = p_raster->getDataset()->ExecuteSQL(queryString.c_str(), nullptr, "SQLITE");

	//OGRPoints to be written, if filename is set
	std::vector<OGRPoint> points;

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

			//only add (and potentially plot/write) a point if it is within the grid polygon
			if (x >= xMin && x <= xMax && y >= yMin && y <= yMax) {
				wktPoints.push_back(point.exportToWkt());
				if (plot) {
					xCoords.push_back(x);
					yCoords.push_back(y);
				}
				if (filename != "") {
					points.push_back(point);
				}
			}

			//set grid vector to be plot
			if (plot) {
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

	if (filename != "") {
		try {
			writeSamplePoints(points, filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}	

	return {wktPoints, {xCoords, yCoords}, grid};
}

PYBIND11_MODULE(systematic, m) {
	m.def("systematic_cpp", &systematic);
}
