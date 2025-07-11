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
	std::mt19937::result_type seed = time(nullptr);
	auto rng = std::bind(
		std::uniform_real_distribution<double>(0, yDiff * xDiff),
		std::mt19937(seed)
	);

	double yCoord = rng() % (yDiff);
	double xCoord = rng() % (xDiff);

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
	std::string extentPolygon = "SELECT 'POLYGON (( " 
		+ std::to_string(xMin) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMax) + ", "
		+ std::to_string(xMax) + " " + std::to_string(yMin) + ", "
		+ std::to_string(xMin) + " " + std::to_string(yMin) + " ))'";

	std::string queryString = "SELECT " + gridFunction + "(ST_GeomFromText("
		+ extentPolygon
		+ "), " + std::to_string(cellSize) + 
		", ST_GeomFromText(" + std::to_string(xCoord) + ", " + std::to_string(yCoord) + "))";

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
		std::cout << p_geometry->getGeometryName() << std::endl;
		for (const auto& p_polygon : *p_geometry->toMultiPolygon()) {
			//generate sample depending on 'location' parameter
			OGRPoint point;
			if (location == "center") {
				p_polygon->(&point);
			}
			else if (location == "corner") {
				point = p_polygon->begin()->begin();
			}
			else { //location == "random"
				OGREnvelope3D *p_envelope;
				p_polygon->getEnvelope(p_envelope);
				double xMin = p_envelope->MinX;
				double xMax = p_envelope->MaxX;
				double xDiff = xMax - xMin;
				double yMin = p_envelope->MinY;
				double yMax = p_envelope->MaxY;
				yDiff = yMax - yMin;

				point.setX(xMin + rng() % xDiff);
				point.setY(yMin + rng() % yDiff);
				while (!p_polygon->Contains(point)) {
					point.setX(xMin + rng() % xDiff);
					point.setY(yMin + rng() % yDiff);
				}
			}

			wktPoints.push_back(point.exportToWkt());

			//if we're plotting, add point to coordinate vectors and plottable vectors to 'grid'
			if (plot) {
				xCoords.push_back(point.getX());
				yCoords.push_back(point.getY());
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
	}

	return {wktPoints, {xCoords, yCoords}, grid};
}

PYBIND11_MODULE(systematic, m) {
	m.def("systematic_cpp", &systematic);
}
