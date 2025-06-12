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
 *
 */
template <typename T, typename U>
std::vector<std::vector<double>> srs_cpp(
	GDALRasterWrapper *p_raster,
	GDALVectorWrapper *p_vector,
	U numSamples)
{
	//Step 1: get dataset
	GDALDataset *p_dataset = p_raster->getDataset();

	//step 2: allocate index array mapping adjusted index to orignial index
	U *p_indexArray = (U *)CPLMalloc(p_raster->getWidth() * p_raster->getHeight());
	U *p_indexArrayUnfilledPointer = p_indexArray;
	U noDataPixelCount = 0;
	
	//step 3: allocate first raster band
	T *p_rasterBand = (T *)CPLMalloc(
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
		(void *)p_rasterBand,			//void *pData
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
	for (U i = 0; i < (U)p_raster->getWidth() * (U)p_raster->getHeight(); i++) {
		T val = p_rasterBand[i];
		if (std::isnan(val) || (double)val == noDataValue) {
			//Step 4.1: increment noDataPixelCount if encountered noData
			noDataPixelCount++;
		}
		else {
			//Step 4.2: add data index to indexArray
			*p_indexArrayUnfilledPointer = i;
			p_indexArrayUnfilledPointer++;
		}
	}
	U numDataPixels = (U)p_raster->getWidth() * (U)p_raster->getHeight() - noDataPixelCount;

	//Step 5: free no longer required raster band
	CPLFree(p_rasterBand);

	//Step 6: generate random number generator using mt19937
	std::mt19937::result_type seed = time(0);
	auto rng = std::bind(
		std::uniform_int_distribution<U>(0, numDataPixels - 1),
		std::mt19937(seed)
	);

	//Step 7: generate numSamples random numbers of data pixels
	std::set<U> samplePixels = {};
	while(samplePixels.size() < numSamples) {
		samplePixels.insert(rng());
	}

	//Step 8: generate coordinate points for each sample index
	double *GT = p_raster->getGeotransform();
	std::vector<double> xCoords;
	std::vector<double> yCoords;
	for( auto samplePixel : samplePixels ) {
		U index = p_indexArray[samplePixel];
		double x_index = index / p_raster->getWidth();
		double y_index = index - (x_index * p_raster->getWidth());
		xCoords.push_back(GT[0] + x_index * GT[1] + y_index * GT[2]);
		yCoords.push_back(GT[3] + x_index * GT[4] + y_index * GT[5]);
	}

	//Step 9: free no longer required index array
	CPLFree(p_indexArray);

	//Step 10: return coordinates
	return {xCoords, yCoords};
}

PYBIND11_MODULE(srs, m) {
	m.def("srs-GDT_Int8-unsigned_short", &srs_cpp<int8_t, unsigned short>);
	m.def("srs-GDT_Int8-unsigned", &srs_cpp<int8_t, unsigned>);
	m.def("srs-GDT_Int8-unsigned_long", &srs_cpp<int8_t, unsigned long>);
	m.def("srs-GDT_Int8-unsigned_long_long", &srs_cpp<int8_t, unsigned long long>);
	
	m.def("srs-GDT_UInt16-unsigned_short", &srs_cpp<uint16_t, unsigned short>);
	m.def("srs-GDT_UInt16-unsigned", &srs_cpp<uint16_t, unsigned>);
	m.def("srs-GDT_UInt16-unsigned_long", &srs_cpp<uint16_t, unsigned long>);
	m.def("srs-GDT_UInt16-unsigned_long_long", &srs_cpp<uint16_t, unsigned long long>);
	
	m.def("srs-GDT_Int16-unsigned_short", &srs_cpp<int16_t, unsigned short>);
	m.def("srs-GDT_Int16-unsigned", &srs_cpp<int16_t, unsigned>);
	m.def("srs-GDT_Int16-unsigned_long", &srs_cpp<int16_t, unsigned long>);
	m.def("srs-GDT_Int16-unsigned_long_long", &srs_cpp<int16_t, unsigned long long>);
	
	m.def("srs-GDT_UInt32-unsigned_short", &srs_cpp<uint32_t, unsigned short>);
	m.def("srs-GDT_UInt32-unsigned", &srs_cpp<uint32_t, unsigned>);
	m.def("srs-GDT_UInt32-unsigned_long", &srs_cpp<uint32_t, unsigned long>);
	m.def("srs-GDT_UInt32-unsigned_long_long", &srs_cpp<uint32_t, unsigned long long>);
	
	m.def("srs-GDT_Int32-unsigned_short", &srs_cpp<int32_t, unsigned short>);
	m.def("srs-GDT_Int32-unsigned", &srs_cpp<int32_t, unsigned>);
	m.def("srs-GDT_Int32-unsigned_long", &srs_cpp<int32_t, unsigned long>);
	m.def("srs-GDT_Int32-unsigned_long_long", &srs_cpp<int32_t, unsigned long long>);
	
	m.def("srs-GDT_Float32-unsigned_short", &srs_cpp<float, unsigned short>);
	m.def("srs-GDT_Float32-unsigned", &srs_cpp<float, unsigned>);
	m.def("srs-GDT_Float32-unsigned_long", &srs_cpp<float, unsigned long>);
	m.def("srs-GDT_Float32-unsigned_long_long", &srs_cpp<float, unsigned long long>);
	
	m.def("srs-GDT_Float64-unsigned_short", &srs_cpp<double, unsigned short>);
	m.def("srs-GDT_Float64-unsigned", &srs_cpp<double, unsigned>);
	m.def("srs-GDT_Float64-unsigned_long", &srs_cpp<double, unsigned long>);
	m.def("srs-GDT_Float64-unsigned_long_long", &srs_cpp<double, unsigned long long>);
}
