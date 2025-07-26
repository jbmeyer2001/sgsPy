/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of random sampling
 * Author: Joseph Meyer
 * Date: June, 2025

 *
 ******************************************************************************/

#include <iostream>
#include <random>

//sgs/utils cpp code
#include "access.h"
#include "raster.h"
#include "vector.h"
#include "write.h"

/**
 * This function uses random sampling to determine the location
 * of sample plots given a raster image which may contain nodata
 * pixels.
 *
 * A single raster band is read from the raster, and each pixel
 * is checked to ensure a sample is not located on a nodata
 * pixel, and the pixel is accessable if access information is provided. 
 * The indeces of the data pixels are saved in a vector. 
 * When all pixels have been read, the indexes are randomly 
 * drawn from the data pixels, converted to geographic coordinates 
 * using the geotransform, and returned as a new GDALVectorWrapper object.
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
 * @returns std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *>
 * 	coordinate and dataset representation of samples
 */
template <typename T, typename U>
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *> 
srs(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string filename)
{
	//Step 1: get dataset and geotransform
	GDALDataset *p_dataset = p_raster->getDataset();
	double *GT = p_raster->getGeotransform();

	//step 2: allocate index array which maps the adjusted index to the orignial index
	std::vector<U> indexes;
	U noDataPixelCount = 0;
	
	//step 3: get first raster band
	T *p_rasterBand = (T *)p_raster->getRasterBand(0);	//GDALRasterWrapper bands are 0 indexed

	//step 4: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
 		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using band
		p_mask = maskInfo.second; //pointer to mask
	}

	//Step 5: iterate through raster band
	double noDataValue = p_dataset->GetRasterBand(1)->GetNoDataValue();
	for (U i = 0; i < (U)p_raster->getWidth() * (U)p_raster->getHeight(); i++) {
		if (p_access) {
			//step 5.1: if the current pixel is not accessable, mark it as nodata and don't read it
			if (((uint8_t *)p_mask)[i] == 0) {
				noDataPixelCount++;
				continue;
			}
		}

		T val = p_rasterBand[i];
		if (std::isnan(val) || (double)val == noDataValue) {
			//Step 5.2: increment noDataPixelCount if encountered noData
			noDataPixelCount++;
		}
		else {
			//Step 5.3: add data index to indexArray
			indexes.push_back(i);
		}
	}
	if (p_access) {
		free(p_accessMaskDataset);
	}

	U numDataPixels = (U)p_raster->getWidth() * (U)p_raster->getHeight() - noDataPixelCount;

	//Step 6: generate random number generator using mt19937	
	std::mt19937::result_type seed = time(nullptr);
	auto rng = std::bind(
		std::uniform_int_distribution<U>(0, numDataPixels - 1),
		std::mt19937(seed)
	);

	//Step 7: generate numSamples random numbers of data pixels, and backup sample pixels if mindist > 0
	//use std::set because we want to iterate in-order because it will be faster
	std::unordered_set<U> samplePixels = {};
	std::unordered_set<U> dontSamplePixels = {};	
	U samplePixelsSize = std::min((mindist == 0.0) ? numSamples : numSamples * 3, (size_t)numDataPixels);

	if (samplePixelsSize > numDataPixels / 2) {
		while (dontSamplePixels.size() < numDataPixels - samplePixelsSize) {
			dontSamplePixels.insert(rng());
		}
		std::shuffle(indexes.begin(), indexes.end(), std::mt19937(seed));
	}
	else {
		while (samplePixels.size() < samplePixelsSize) {
			samplePixels.insert(rng());
		}
	}
	
	//step 8: create new in-memory dataset to store sample points
	//TODO error check this?
	GDALAllRegister();
	GDALDataset *p_sampleDataset = GetGDALDriverManager()->GetDriverByName("MEM")->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	OGRLayer *p_layer = p_sampleDataset->CreateLayer("samples", nullptr, wkbPoint, nullptr);

	//Step 9: generate coordinate points for each sample index, and only add if they're outside of mindist
	std::vector<double> xCoords;
	std::vector<double> yCoords;
	size_t pointsAdded = 0;

	if (dontSamplePixels.size() == 0) {
		for( auto samplePixel : samplePixels ) {
			U index = indexes[samplePixel];	
			double yIndex = index / p_raster->getWidth();
			double xIndex = index - (yIndex * p_raster->getWidth());
			double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
			double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
			OGRPoint newPoint = OGRPoint(xCoord, yCoord);
		
			if (mindist != 0.0 && pointsAdded != 0) {
				bool add = true;
				for (const auto &p_feature : *p_layer) {
					OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
					if (newPoint.Distance(p_point) < mindist) {
						add = false;
						break;
					}
				}
				if (!add) {
					continue;
				}
			}
			
			OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
			p_feature->SetGeometry(&newPoint);
			p_layer->CreateFeature(p_feature);
			OGRFeature::DestroyFeature(p_feature);

			pointsAdded++;
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
	
			if (pointsAdded == numSamples) {
				break;
			}
		}
	}
	else {
		for (U i = 0; i < indexes.size(); i++) {
			if (dontSamplePixels.find(i) != dontSamplePixels.end()) {
				continue;
			}	

			U index = indexes[i];
			double yIndex = index / p_raster->getWidth();
			double xIndex = index - (yIndex * p_raster->getWidth());
			double yCoord = GT[3] + xIndex * GT[4] + yIndex * GT[5];
			double xCoord = GT[0] + xIndex * GT[1] + yIndex * GT[2];
			OGRPoint newPoint = OGRPoint(xCoord, yCoord);
		
			if (mindist != 0.0 && pointsAdded != 0) {
				bool add = true;
				for (const auto &p_feature : *p_layer) {
					OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
					if (newPoint.Distance(p_point) < mindist) {
						add = false;
						break;
					}
				}
				if (!add) {
					continue;
				}
			}
				
			OGRFeature *p_feature = OGRFeature::CreateFeature(p_layer->GetLayerDefn());
			p_feature->SetGeometry(&newPoint);
			p_layer->CreateFeature(p_feature);
			OGRFeature::DestroyFeature(p_feature);
			
			pointsAdded++;
			xCoords.push_back(xCoord);
			yCoords.push_back(yCoord);
		
			if (pointsAdded == numSamples) {
				break;
			}		
		}
	}

	//step 10: create GDALVectorWrapper with dataset containing points
	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_sampleDataset);

	//Step 11: write vector of points if given filename
	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	return {{xCoords, yCoords}, p_sampleVectorWrapper};
}

/**
 * Having template types which rely on dynamic information (such as
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
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *> 
srsTypeSpecifier(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string filename) 
{
	std::string minIndexIntType = p_raster->getMinIndexIntType(true); //singleBand = true 
	switch (p_raster->getRasterType()) {
		case GDT_Int8:
		if(minIndexIntType == "unsigned_short") { return srs<int8_t, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<int8_t, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<int8_t, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<int8_t, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		case GDT_UInt16:
		if(minIndexIntType == "unsigned_short") { return srs<uint16_t, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<uint16_t, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<uint16_t, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<uint16_t, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		break;
		case GDT_Int16:
		if(minIndexIntType == "unsigned_short") { return srs<int16_t, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<int16_t, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<int16_t, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<int16_t, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		break;
		case GDT_UInt32:
		if(minIndexIntType == "unsigned_short") { return srs<uint32_t, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<uint32_t, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<uint32_t, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<uint32_t, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		break;
		case GDT_Int32:
		if(minIndexIntType == "unsigned_short") { return srs<int32_t, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<int32_t, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<int32_t, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<int32_t, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		break;
		case GDT_Float32:
		if(minIndexIntType == "unsigned_short") { return srs<float, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<float, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<float, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<float, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		break;
		case GDT_Float64:
		if(minIndexIntType == "unsigned_short") { return srs<double, unsigned short>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned") { return srs<double, unsigned>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long") { return srs<double, unsigned long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		if(minIndexIntType == "unsigned_long_long") { return srs<double, unsigned long long>(p_raster, numSamples, mindist, p_access, layerName, buffInner, buffOuter, filename); }
		break;
		default: 
		throw std::runtime_error("GDALDataType not one of the accepted types.");
	}
	throw std::runtime_error("type " + minIndexIntType + " not a valid type.");
}

//TODO: maybe find a better (less verbose) way of optional Python arguments than this...

/*
 * srs function for when access has not been specified. 
 * Calls srsTypeSpecifier() which calls srs().
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *> 
srs_cpp(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	double mindist,
	std::string filename) 
{
	return srsTypeSpecifier(
		p_raster, 
		numSamples,
		mindist, 
		nullptr, 
		"",
		0, 
		0, 
		filename
	);
}

/*
 * srs function for when access has been specified.
 * Calls srsTypeSpecifier() which calls srs().
 */
std::pair<std::vector<std::vector<double>>, GDALVectorWrapper *> 
srs_cpp_access(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	std::string filename) 
{
	return srsTypeSpecifier(
		p_raster, 
		numSamples,
		mindist, 
		p_access,
	       	layerName,	
		buffInner, 
		buffOuter, 
		filename
	);
}

PYBIND11_MODULE(srs, m) {
	m.def("srs_cpp", &srs_cpp);
	m.def("srs_cpp_access", &srs_cpp_access);
}
