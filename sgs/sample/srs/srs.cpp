/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of random sampling
 * Author: Joseph Meyer
 * Date: June, 2025

 *
 ******************************************************************************/

#include <random>

//sgs/utils cpp code
#include "raster.h"
#include "vector.h"

/**
 * This function uses random sampling to determine the location
 * of sample plots given a raster image which may contain nodata
 * pixels.
 *
 * //TODO add access vector information when functionality
 * is implemented.
 *
 * A single raster band is read from the raster, and each pixel
 * is checked to ensure a sample is not located on a nodata
 * pixel. The indeces of the data pixels are saved in another 
 * allocated matrix. When all pixels have been read, the indexes
 * are randomly drawn from the data pixels, converted to 
 * geographic coordinates using the geotransform, and returned.
 *
 * The function is a template function which contains the data
 * type of the raster (T), as well as the type of the 
 * index matrix (U). The use of U as a template argument is
 * to ensure the index matrix data type can represent every
 * index in the raster without overflow, and do it with
 * the most efficient use of memory.
 *
 * @param GDALRasterWrapper * a pointer to the raster image we're sampling
 * @param GDALRasterVector * a pointer to an access vector image if it exists
 * @param U <unsigned short, unsigned, unsigned long, unsigned long long> the number of samples
 * @returns std::pair<std::vector<std::vector<double>>, std::vector<std::string>>
 * 		coordinate and wkt representation of samples
 */
template <typename T, typename U>
std::pair<std::vector<std::vector<double>>, std::vector<std::string>> 
srs(
	GDALRasterWrapper *p_raster,
	double mindist,
	//GDALVectorWrapper *p_vector,
	U numSamples,
	std::string filename)
{
	//Step 1: get dataset
	GDALDataset *p_dataset = p_raster->getDataset();

	//step 2: allocate index array mapping adjusted index to orignial index
	U *p_indexArray = (U *)std::malloc(p_raster->getWidth() * p_raster->getHeight() * sizeof(U));
	U *p_indexArrayUnfilledPointer = p_indexArray;
	U noDataPixelCount = 0;
	
	//step 3: allocate first raster band
	T *p_rasterBand = (T *)std::malloc(
		p_raster->getWidth() * 		//width
		p_raster->getHeight() * 	//height
		p_raster->getRasterTypeSize()	//num bytes per pixel
	);
	CPLErr err = p_dataset->GetRasterBand(1)->RasterIO(
		GF_Read,			//GDALRWFlag eRWFlag
		0,				//int nXOff
		0,				//int nYOff
		p_raster->getWidth(),		//int nXSize
		p_raster->getHeight(),		//int nYSize
		(void *)p_rasterBand,		//void *pData
		p_raster->getWidth(),		//int nBufXSize
		p_raster->getHeight(),		//int nBufYSize
		p_raster->getRasterType(),	//GDALDataType eBufType
		0,				//int nPixelSpace
		0				//int nLineSpace
	);
	if (err) {
		throw std::runtime_error("Error reading dataset band.");
	}

	//Step 4: iterate through raster band
	double noDataValue = p_dataset->GetRasterBand(1)->GetNoDataValue();
	U j = 0;
	for (U i = 0; i < (U)p_raster->getWidth() * (U)p_raster->getHeight(); i++) {
		T val = p_rasterBand[i];
		if (std::isnan(val) || (double)val == noDataValue) {
			//Step 4.1: increment noDataPixelCount if encountered noData
			noDataPixelCount++;
		}
		else {
			//Step 4.2: add data index to indexArray
			p_indexArrayUnfilledPointer[j] = i;
			j++;
		}
	}
	U numDataPixels = (U)p_raster->getWidth() * (U)p_raster->getHeight() - noDataPixelCount;

	//step 5: free no longer required raster band
	CPLFree(p_rasterBand);

	//Step 6: generate random number generator using mt19937	
	std::mt19937::result_type seed = time(nullptr);
	auto rng = std::bind(
		std::uniform_int_distribution<U>(0, numDataPixels - 1),
		std::mt19937(seed)
	);

	//Step 7: generate numSamples random numbers of data pixels, and backup sample pixels if mindist > 0
	//use std::set because we want to iterate in-order because it will be faster
	std::set<U> samplePixels = {}; 	
	
	//use std::unordered_set because we want to iterate out-of-order because it is required by the implementation
	std::unordered_set<U> backupSamplePixels = {};

	while(samplePixels.size() < numSamples) {
		samplePixels.insert(rng());
	}
	if (mindist != 0.0)  {
		while(backupSamplePixels.size() < numSamples) {
			U pixel = rng();
			if (samplePixels.find(pixel) == samplePixels.end()) {
				backupSamplePixels.insert(pixel);
			}
		}
	}

	//Step 8: generate coordinate points for each sample index, and only add if they're outside of mindist
	double *GT = p_raster->getGeotransform();
	std::vector<double> xCoords;
	std::vector<double> yCoords;
	std::vector<OGRPoint> points;
	std::vector<std::string> wktPoints;

	for( auto samplePixel : samplePixels ) {
		U index = p_indexArray[samplePixel];
		double xIndex = index / p_raster->getWidth();
		double yIndex = index - (xIndex * p_raster->getWidth());
		double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
		double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
		OGRPoint newPoint = OGRPoint(xCoord, yCoord);
		
		if (mindist != 0.0 || points.size() == 0) {
			U pIndex = 0;
			while ((newPoint.Distance(&points[pIndex]) > mindist) && (pIndex < points.size())) {
				pIndex++;
			}
			if (pIndex != points.size()) {
				continue;
			}
		}

		points.push_back(newPoint);
		wktPoints.push_back(newPoint.exportToWkt());
		xCoords.push_back(xCoord);
		yCoords.push_back(yCoord);
	}

	//if we need more sample points, iterate through the backups until we have
	//either run out of backups or have the required number of samples
	if (mindist != 0.0 && points.size() < numSamples) {
		for ( auto samplePixel : backupSamplePixels ) {
			U index = p_indexArray[samplePixel];
			double xIndex = index / p_raster->getWidth();
			double yIndex = index - (xIndex * p_raster->getWidth());
			double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
			double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
			OGRPoint newPoint = OGRPoint(xCoord, yCoord);

			U pIndex = 0;
			while ((newPoint.Distance(&points[pIndex]) > mindist) && (pIndex < points.size())) {
				pIndex++;
			}

			if (pIndex == points.size()) {
				points.push_back(newPoint);
				wktPoints.push_back(newPoint.exportToWkt());
				xCoords.push_back(xCoord);
				yCoords.push_back(yCoord);
			}

			if (points.size() == numSamples) { 
				break;
			}
		}
	}

	//Step 9: free no longer required index array
	CPLFree(p_indexArray);

	//Step 10: return as coordinates and wkt
	return {{xCoords, yCoords}, wktPoints};
}

/**
 * Having template types which rely in dynamic information (such as
 * the pixel type of the added raster, or the number of pixels in 
 * the raster) require an unfortunate amount of boilerplate code.
 *
 * This is an attempt to condense as much of the annoying boilerplate
 * into a single place.
 *
 * This function uses type information of the raster pixel type,
 * as well as the minimally sized unsigned int type which can represent
 * all necessary indices.
 *
 * A call is made to srs() with the necessary data type template
 * arguments depending on raster parameters.
 *
 * @returns std::pair<std::vector<std::vector<double>>, std::vector<std::string>> 
 * 		coordinate and wkt representation of samples
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::string>> 
srsTypeSpecifier(
	GDALRasterWrapper *p_raster,
	double mindist,
	//GDALVectorWrapper *p_vector,
	size_t numSamples,
	std::string filename) 
{
	std::string minIndexIntType = p_raster->getMinIndexIntType(true); //singleBand = true 
	switch (p_raster->getRasterType()) {
		case GDT_Int8:
		if(minIndexIntType == "unsigned_short") { return srs<int8_t, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<int8_t, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<int8_t, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<int8_t, unsigned long long>(p_raster, mindist, numSamples, filename); }
		case GDT_UInt16:
		if(minIndexIntType == "unsigned_short") { return srs<uint16_t, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<uint16_t, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<uint16_t, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<uint16_t, unsigned long long>(p_raster, mindist, numSamples, filename); }
		break;
		case GDT_Int16:
		if(minIndexIntType == "unsigned_short") { return srs<int16_t, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<int16_t, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<int16_t, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<int16_t, unsigned long long>(p_raster, mindist, numSamples, filename); }
		break;
		case GDT_UInt32:
		if(minIndexIntType == "unsigned_short") { return srs<uint32_t, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<uint32_t, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<uint32_t, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<uint32_t, unsigned long long>(p_raster, mindist, numSamples, filename); }
		break;
		case GDT_Int32:
		if(minIndexIntType == "unsigned_short") { return srs<int32_t, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<int32_t, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<int32_t, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<int32_t, unsigned long long>(p_raster, mindist, numSamples, filename); }
		break;
		case GDT_Float32:
		if(minIndexIntType == "unsigned_short") { return srs<float, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<float, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<float, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<float, unsigned long long>(p_raster, mindist, numSamples, filename); }
		break;
		case GDT_Float64:
		if(minIndexIntType == "unsigned_short") { return srs<double, unsigned short>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned") { return srs<double, unsigned>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<double, unsigned long>(p_raster, mindist, numSamples, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<double, unsigned long long>(p_raster, mindist, numSamples, filename); }
		break;
		default: 
		throw std::runtime_error("GDALDataType not one of the accepted types.");
	}
	throw std::runtime_error("type " + minIndexIntType + " not a valid type.");
}

PYBIND11_MODULE(srs, m) {
	m.def("srs_cpp", &srsTypeSpecifier);
}
