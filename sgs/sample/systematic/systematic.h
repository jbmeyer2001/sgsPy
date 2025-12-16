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

#include "utils/access.h"
#include "utils/existing.h"
#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

namespace systematic {

/**
 * Helper function for generating a vector geometry containing polygons of the accessible
 * area. This function is called if an access vector is given.
 *
 * The input vector dataset must be comprised of solely LineString and MultiLineString geometries.
 * The LineStrings are buffered by buff outer, then have buff_inner subtracted from them.
 *
 * The returned geometry is checked to see whether it contains points which are trying to be
 * sampled. If the geometry does not contain a given point, that point is considered inaccessible
 * and thus not sampled.
 *
 * @param GDALVectorWrapper *p_access
 * @param std::string layerName
 * @param double buffInner
 * @param double buffOuter
 */
OGRGeometry *getAccessPolygon(GDALVectorWrapper *p_access, std::string layerName, double buffInner, double buffOuter) {
	//step 1: create multipolygon buffers
	OGRMultiPolygon *buffInnerPolygons = new OGRMultiPolygon;
	OGRMultiPolygon *buffOuterPolygons = new OGRMultiPolygon;
	OGRGeometry *p_polygonMask;

	//step 2: add geometries to access polygon buffers
	for (const auto& p_feature : *p_access->getLayer(layerName.c_str())) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		OGRwkbGeometryType type = wkbFlatten(p_geometry->getGeometryType());

		switch (type) {
			case OGRwkbGeometryType::wkbLineString: {
				buffOuterPolygons->addGeometry(p_geometry->Buffer(buffOuter));
				if (buffInner != 0) {
					buffInnerPolygons->addGeometry(p_geometry->Buffer(buffInner));
				}
				break;
			}
			case OGRwkbGeometryType::wkbMultiLineString: {
				for (const auto& p_lineString : *p_geometry->toMultiLineString()) {
					buffOuterPolygons->addGeometry(p_lineString->Buffer(buffOuter));
					if (buffInner != 0) {
						buffInnerPolygons->addGeometry(p_lineString->Buffer(buffInner));
					}
				}
				break;
			}
			default: {
				throw std::runtime_error("access polygon geometry type must be LineString or MultiLineString");
			}
		}
	}	

	//step 3: generate the polygon mask and free no longer used memory
	if (buffInner == 0) {
		p_polygonMask = buffOuterPolygons->UnionCascaded();
		free(buffOuterPolygons);
	}
	else {
		OGRGeometry *buffOuterUnion = buffOuterPolygons->UnionCascaded();
		OGRGeometry *buffInnerUnion = buffInnerPolygons->UnionCascaded();
		p_polygonMask = buffOuterUnion->Difference(buffInnerUnion);
		free(buffOuterPolygons);
		free(buffInnerPolygons);
		free(buffOuterUnion);
		free(buffInnerUnion);
	}

	return p_polygonMask;
}

/**
 * Helper function for ensuring the x and y coordinates are within the raster extent.
 *
 * @param double x
 * @param double y
 * @param double xMin
 * @param double xMax
 * @param double yMin
 * @param double yMax
 */
inline bool
checkExtent(double x, double y, double xMin, double xMax, double yMin, double yMax) {
	return (x >= xMin && x <= xMax && y >= yMin && y <= yMax); 
}

/**
 * Helper function for checking to see whether a pixel occurs within accessible area.
 *
 * @param OGRPoint *p_point
 * @param OGRGeometry *p_geometry
 */
inline bool
checkAccess(OGRPoint *p_point, OGRGeometry *p_geometry) {
	return !p_geometry || p_point->Within(p_geometry);
}

/**
 * Helper function for checking to see whether a pixel is already an existing sample location.
 *
 * @param double x
 * @param double y
 * @param Existing& existing
 */
inline bool
checkExisting(double x, double y, Existing& existing) {
	return !existing.used || !existing.containsCoordinates(x, y);
}

/**
 * Helper function for checking to see whether a coordinate occurs in an area of nodata.
 *
 * the 'force' parameter is first checked, because if force is not true, then samples
 * are allowed to occur on nodata pixels.
 *
 * THe inverted geotransform is used to find the x and y values from the coordinates
 * given. Then, every raster band within the input raster is checked, and if that 
 * pixel is a no data value in any of them false is returned. the GDALRasterBand
 * RasterIO function is used to read the desried pixel from the raster.
 *
 * @param GDALRasterWrapper *p_raster,
 * @param double *IGT
 * @param double xCoord
 * @param double yCoord
 * @param bool force
 */
inline bool
checkNotNan(GDALRasterWrapper *p_raster, double *IGT, double xCoord, double yCoord, bool force) {
	if (force) {
		int x = static_cast<int>(IGT[0] + xCoord * IGT[1] + yCoord * IGT[2]);
		int y = static_cast<int>(IGT[3] + xCoord * IGT[4] + yCoord * IGT[5]);
	
		for (int i = 0; i < p_raster->getBandCount(); i++) {
			GDALRasterBand *p_band = p_raster->getRasterBand(i);
			
			double val;
			p_band->RasterIO(GF_Read, x, y, 1, 1, &val, 1, 1, GDT_Float64, 0, 0);
	
			if (val == p_band->GetNoDataValue() || std::isnan(val)) {
				return false;
			}	
		}	
	}

	return true;
}	

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
 * location parameter (centers, corners, or random), and the samples are
 * saved to an OGRLayer which is part of GDALDataset. Plot-required data
 * is saved if plot is true (to later be utilized by the Python side
 * of the application with matplotlib).
 *
 * If the an access vector is given, then polygons of the accessible area
 * are created, and each sample is check to ensure it falls within the
 * accessible area.
 *
 * If an existing vector is given, all of the sample points within the
 * existing vector are added, and each point is checked to ensure it
 * has not already been added by virtue of already existing as a sample
 * point.
 *
 * If the force parameter is given, every sample added is checked against
 * the input raster to ensure the sample does not fall in a no data pixel.
 * If it does the sample is thrown out. In the case where each grid cell
 * is randomly sampled, 10 tries are allowed to find a point which is
 * contains a data pixel otherwise that cell is not sampled.
 *
 * @param GDALRasterWrapper *p_raster raster to be systematically sampled
 * @param double cellSize the size of the grid cell shapes
 * @param std::string shape the shape of the grid cells
 * @param std::string location the location within a cell to sample
 * @param bool plot whether to save and return plot-required data
 * @param std::string filename to write to or "" if not to write
 * @returns std::tuple<
 * 		GDALVectorWrapper *,
 * 		std::vector<std::vector<double>>,
 * 		std::vector<std::vector<std::vector<double>>>
 * 	>
 * 	Wrapper containing GDALDataset vector of sample points,
 * 	the 2d vector of doubles contains sample points to plot
 * 	the 3d vector of doubles contains grid cells to plot
 */
std::tuple<
	GDALVectorWrapper *, //GDALDataset containing sample points
	std::vector<std::vector<double>>, //array of samples to plot
	std::vector<std::vector<std::vector<double>>> //array of grid to plot
>
systematic(
	GDALRasterWrapper *p_raster,
	double cellSize,
	std::string shape,
	std::string location,
	GDALVectorWrapper *p_existing,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool force,
	bool plot,
	std::string filename)
{
	GDALAllRegister();
	
	double *GT; 	//geotransform
	double IGT[6];	//inverse geotransform
	GT = p_raster->getGeotransform();
	GDALInvGeoTransform(GT, IGT);	

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

	//query to create grid
	OGRLayer *p_gridTest = p_raster->getDataset()->ExecuteSQL(queryString.c_str(), nullptr, "SQLITE");

	//create output dataset before anything that might take a long time (in case creation of dataset fails)
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
	if(!p_driver) {
		throw std::runtime_error("unable to create output sample dataset driver.");
	}
	GDALDataset *p_sampleDataset = p_driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	if (!p_sampleDataset) {
		throw std::runtime_error("unable to create output dataset with driver.");
	}

	GDALVectorWrapper *p_wrapper = new GDALVectorWrapper(p_sampleDataset, std::string(p_raster->getDataset()->GetProjectionRef()));
	OGRLayer *p_sampleLayer = p_sampleDataset->CreateLayer("samples", p_wrapper->getSRS(), wkbPoint, nullptr);	
	if (!p_sampleLayer) {
		throw std::runtime_error("unable to create output dataset layer.");
	}

	//get access polygon if access is given
	OGRGeometry *access = nullptr;
	if (p_access) {
		//occasionally, samples end up being placed just outside the accessible area,
		//typically when the buffer sizes are multiples of the pixel size. 
		//
		//Adjusting the buffers minorly fixes this problem.
		double pixelSize = std::min(p_raster->getPixelHeight(), p_raster->getPixelWidth());
		buffOuter = buffOuter - pixelSize / 50;
		buffInner = buffInner == 0 ? 0 : buffInner + pixelSize / 50;

		access = getAccessPolygon(p_access, layerName, buffInner, buffOuter);
	}

	//coordinate representations the samples as vectors only if PLOT is true, returned to the (Python) caller
	std::vector<double> xCoords, yCoords;

	//create existing struct
	Existing existing(p_existing, GT, p_raster->getWidth(), p_sampleLayer, plot, xCoords, yCoords);

	//grid represents the grid used to create the sampels only if PLOT is true, and is returned to the (Python) caller
	std::vector<std::vector<std::vector<double>>> grid;

	//iterate through the polygons in the grid and populate the ruturn data depending on user inputs
	for (const auto& p_feature : *p_gridTest) {
		OGRGeometry *p_geometry = p_feature->GetGeometryRef();
		for (const auto& p_polygon : *p_geometry->toMultiPolygon()) {
			//generate sample depending on 'location' parameter
			OGRPoint point;
			OGRPoint secondPoint;
			if (location == "centers") {
				p_polygon->Centroid(&point);

				double x = point.getX();
				double y = point.getY();
				if (checkExtent(x, y, xMin, xMax, yMin, yMax) &&
				    checkAccess(&point, access) &&
				    checkExisting(x, y, existing) &&
				    checkNotNan(p_raster, IGT, x, y, force)) 
				{
					addPoint(&point, p_sampleLayer);

					if (plot) {
						xCoords.push_back(x);
						yCoords.push_back(y);
					}	
				}
			}
			else if (location == "corners") {
				(*p_polygon->begin())->getPoint(0, &point);
				(*p_polygon->begin())->getPoint(1, &secondPoint);

				double x = point.getX();
				double y = point.getY();
				if (checkExtent(x, y, xMin, xMax, yMin, yMax) &&
				    checkAccess(&point, access) &&
				    checkExisting(x, y, existing) &&
				    checkNotNan(p_raster, IGT, x, y, force)) 
				{
					addPoint(&point, p_sampleLayer);

					if (plot) {
						xCoords.push_back(x);
						yCoords.push_back(y);
					}	
				}

				x = secondPoint.getX();
				y = secondPoint.getY();
				if (checkExtent(x, y, xMin, xMax, yMin, yMax) &&
				    checkAccess(&secondPoint, access) &&
				    checkExisting(x, y, existing) &&
				    checkNotNan(p_raster, IGT, x, y, force)) 
				{
					addPoint(&secondPoint, p_sampleLayer);

					if (plot) {
						xCoords.push_back(x);
						yCoords.push_back(y);
					}	
				}

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

				int tries = 0;
				bool found = false;

				while (tries < 10 && !found) {
					point.setX(xMinEnv + rng() / (rngMax / xDiffEnv));
					point.setY(yMinEnv + rng() / (rngMax / yDiffEnv));
					while (!p_polygon->Contains(&point)) {
						point.setX(xMinEnv + rng() / (rngMax / xDiffEnv));
						point.setY(yMinEnv + rng() / (rngMax / yDiffEnv));
					}

					double x = point.getX();
					double y = point.getY();
					if (checkExtent(x, y, xMin, xMax, yMin, yMax) &&
				    	    checkAccess(&point, access) &&
				    	    checkExisting(x, y, existing) &&
					    checkNotNan(p_raster, IGT, x, y, force)) 
					{
						found = true;
						addPoint(&point, p_sampleLayer);

						if (plot) {
							xCoords.push_back(x);
							yCoords.push_back(y);
						}	
					}

					tries++;
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
			p_wrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}	

	return {p_wrapper, {xCoords, yCoords}, grid};
}

} //namespace systematic
