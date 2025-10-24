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
std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t> 
srs(
	GDALRasterWrapper *p_raster,
	size_t numSamples,
	double mindist,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string filename)
{
	//Step 1: get dataset and geotransform
	GDALDataset *p_dataset = p_raster->getDataset();
	double *GT = p_raster->getGeotransform();

	//step 2: allocate index array which maps the adjusted index to the orignial index
	std::vector<size_t> indexes;
	size_t noDataPixelCount = 0;
	
	//step 3: get first raster band
	GDALRasterBand *p_band = p_raster->getRasterBand(0);
	double noDataValue = p_band->GetNoDataValue();
	GDALDataType type = p_raster->getRasterBandType(0);
	void *p_data = VSIMalloc3(
		p_raster->getHeight(),
		p_raster->getWidth(),
		p_raster->getRasterBandTypeSize(0) //per pixel size of band 1
	);
	CPLErr err = p_band->RasterIO(
		GF_Read,		//GDALRWFlat eRWFlag
		0,			//int nXOff
		0,			//int nYOff
		p_raster->getWidth(),	//int nXSize
		p_raster->getHeight(),	//int nYSize
		p_data,			//void *pData
		p_raster->getWidth(),	//int nBufXSize
		p_raster->getHeight(),	//int nBufYSize
		type,			//GDALDataType eBufType
		0,			//int nPixelSpace
		0			//int nLineSpace
	);
	if (err) {
		throw std::runtime_error("error reading raster band from dataset.");
	}

	//step 4: get access mask if access is defined
	GDALDataset *p_accessMaskDataset = nullptr;
	void *p_mask = nullptr;
	if (p_access) {
		std::pair<GDALDataset *, void *> maskInfo = getAccessMask(p_access, p_raster, layerName, buffInner, buffOuter);
 		p_accessMaskDataset = maskInfo.first; //GDALDataset to free after using band
		p_mask = maskInfo.second; //pointer to mask
	}

	//Step 5: iterate through raster band
	for (size_t i = 0; i < (size_t)p_raster->getWidth() * (size_t)p_raster->getHeight(); i++) {
		if (p_access && ((uint8_t *)p_mask)[i] == 0) {
			//step 5.1: if the current pixel is not accessable, mark it as nodata and don't read it
			noDataPixelCount++;
			continue;
		}

		bool isNan = false;
		switch(type) {
			case GDT_Int8: {
				int8_t val = ((int8_t *)p_data)[i];
				isNan = std::isnan(val) || (double)val == noDataValue;
				break;
			}
			case GDT_UInt16: {
				uint16_t val = ((uint16_t *)p_data)[i];
				isNan = std::isnan(val) || (double)val == noDataValue;
				break;
			}
			case GDT_Int16: {
				int16_t val = ((int16_t *)p_data)[i];
				isNan = std::isnan(val) || (double)val == noDataValue;
				break;
			}
			case GDT_UInt32: {
				uint32_t val = ((uint32_t *)p_data)[i];
				isNan = std::isnan(val) || (double)val == noDataValue;
				break;
			}
			case GDT_Int32: {
				int32_t val = ((int32_t *)p_data)[i];
				isNan = std::isnan(val) || (double)val == noDataValue;
				break;
			}
			case GDT_Float32: {
				float val = ((float *)p_data)[i];
				isNan = std::isnan(val) || (double)val == noDataValue;
				break;
			}
			case GDT_Float64: {
				double val = ((double *)p_data)[i];
				isNan = std::isnan(val) || val == noDataValue;
				break;
			}
			default:
				throw std::runtime_error("raster pixel data type not supported.");
		}

		if (isNan) {
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

	size_t numDataPixels = (size_t)p_raster->getWidth() * (size_t)p_raster->getHeight() - noDataPixelCount;

	//Step 6: generate random number generator using mt19937	
	std::mt19937::result_type seed = time(nullptr);
	auto rng = std::bind(
		std::uniform_int_distribution<size_t>(0, numDataPixels - 1),
		std::mt19937(seed)
	);

	//Step 7: generate numSamples random numbers of data pixels, and backup sample pixels if mindist > 0
	//use std::set because we want to iterate in-order because it will be faster
	std::unordered_set<size_t> samplePixels = {};
	std::unordered_set<size_t> dontSamplePixels = {};	
	size_t samplePixelsSize = std::min((mindist == 0.0) ? numSamples : numSamples * 3, (size_t)numDataPixels);

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
	std::vector<double> xCoords, yCoords;
	size_t pointsAdded = 0;

	if (dontSamplePixels.size() == 0) {
		for( auto samplePixel : samplePixels ) {
			size_t index = indexes[samplePixel];	
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
			
			if (plot) {
				xCoords.push_back(xCoord);
				yCoords.push_back(yCoord);
			}
	
			if (pointsAdded == numSamples) {
				break;
			}
		}
	}
	else {
		for (size_t i = 0; i < indexes.size(); i++) {
			if (dontSamplePixels.find(i) != dontSamplePixels.end()) {
				continue;
			}	

			size_t index = indexes[i];
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

			if (plot) {
				xCoords.push_back(xCoord);
				yCoords.push_back(yCoord);
			}

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

	return {{xCoords, yCoords}, p_sampleVectorWrapper, pointsAdded};
}

PYBIND11_MODULE(srs, m) {
	m.def("srs_cpp", &srs, 
		pybind11::arg("p_raster"),
		pybind11::arg("numSamples"),
		pybind11::arg("mindist"),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("filename"));

}
