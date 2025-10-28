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

#include <xoshiro.h>

#include "access.h"
#include "existing.h"
#include "helper.h"
#include "raster.h"
#include "vector.h"

template <typename T>
inline void
processBlock(
	RasterBandMetaData& band,
	Access& access,
	Existing& existing,
	std::vector<Index>& indices,
	std::vector<bool>& randVals,
	int& randValIndex,
	int xBlock,
	int yBlock,
	int xValid,
	int yValid,
	int64_t& accessible,
	int64_t& notAccessible)
{
	T nan = static_cast<T>(band.nan);
	int8_t *p_access = reinterpret_cast<int8_t *>(access.band.p_buffer);

	for (int y = 0; y < yValid; y++) {
		size_t blockIndex = static_cast<size_t>(y * band.xBlockSize);
		for (int x = 0; x < xValid; x++) {
			//GET VAL
			T val = getPixelValueDependingOnType<T>(band.type, band.p_buffer, blockIndex);

			//CHECK NAN
			bool isNan = val == nan || std::isnan(val);
			if (isNan) {
				blockIndex++;
				continue;
			}

			//CHECK ACCESS
			if (access.used && p_access[blockIndex] != 1) {
				notAccessible++;
				blockIndex++;
				continue;
			}
			accessible++;

			Index index = {x + xBlock * band.xBlockSize, y + yBlock * band.yBlockSize};
			//CHECK EXISTING
			if (existing.used && existing.containsIndex(index.x, index.y)) {
				blockIndex++;
				continue;
			}

			//CHECK RNG
			if (!randVals[randValIndex]) {
				randValIndex++;
				blockIndex++;
				continue;
			}
			randValIndex++;

			//ADD TO INDICES
			indices.push_back(index);

			//INCREMENT WITHIN-BLOCK INDEX
			blockIndex++;
		}
	}
}

/**
 *
 */
inline uint64_t
getProbabilityMultiplier(GDALRasterWrapper *p_raster, int numSamples, bool useMindist, double accessibleArea) {
	double height = static_cast<double>(p_raster->getHeight());
	double width = static_cast<double>(p_raster->getWidth());
	double samples = static_cast<double>(numSamples);

	double numer = samples * 4 * (useMindist ? 3 : 1);
	double denom = height * width;

	if (accessibleArea != -1) {
		double pixelHeight = static_cast<double>(p_raster->getPixelHeight());
		double pixelWidth = static_cast<double>(p_raster->getPixelWidth());
		double totalArea = width * pixelWidth * height * pixelHeight;

		numer *= (totalArea / accessibleArea);
	}

	uint8_t bits = static_cast<uint8_t>(std::ceil(std::log2(denom) - std::log2(numer)));
	return (bits <= 0) ? 0 : (1 << bits) - 1;
}

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
	GDALVectorWrapper *p_existing,
	GDALVectorWrapper *p_access,
	std::string layerName,
	double buffInner,
	double buffOuter,
	bool plot,
	std::string tempFolder,
	std::string filename)
{
	GDALAllRegister();

	int width = p_raster->getWidth();
	int height = p_raster->getHeight();
	double *GT = p_raster->getGeotransform();
	bool useMindist = mindist != 0;
	std::mutex mutex;

	std::vector<double> xCoords, yCoords;
	std::vector<size_t> indexes;

	//step 3: get first raster band
	RasterBandMetaData band;
	band.p_band = p_raster->getRasterBand(0);
	band.type = p_raster->getRasterBandType(0);
	band.size = p_raster->getRasterBandTypeSize(0);
	band.nan = band.p_band->GetNoDataValue();
	band.p_mutex = &mutex;
	band.p_band->GetBlockSize(&band.xBlockSize, &band.yBlockSize);
	band.p_buffer = VSIMalloc3(band.xBlockSize, band.yBlockSize, band.size);

	//create output dataset before doing anything which will take a long time in case of failure.
	GDALDriver *p_driver = GetGDALDriverManager()->GetDriverByName("MEM");
	if (!p_driver) {
		throw std::runtime_error("unable to create output sample dataset driver.");
	}
	GDALDataset *p_samples = p_driver->Create("", 0, 0, 0, GDT_Unknown, nullptr);
	if (!p_samples) {
		throw std::runtime_error("unable to create output dataset with driver.");
	}
	OGRLayer *p_layer = p_samples->CreateLayer("samples", nullptr, wkbPoint, nullptr);
	if (!p_layer) {
		throw std::runtime_error("unable to create output dataset layer.");
	}

	//generate access structure
	Access access(
		p_access,
		p_raster,
		layerName,
		buffInner,
		buffOuter,
		true,
		tempFolder,
		band.xBlockSize,
		band.yBlockSize
	);

	if (access.used) {
		access.band.p_buffer = VSIMalloc3(band.xBlockSize, band.yBlockSize, access.band.size);
	}

	//generate existing structure
	Existing existing(
		p_existing,
		GT,
		p_raster->getWidth(),
		p_layer,
		plot,
		xCoords,
		yCoords
	);

	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (height + band.yBlockSize - 1) / band.yBlockSize;

	std::vector<bool> randVals(band.xBlockSize * band.yBlockSize);
	int randValIndex = band.xBlockSize * band.yBlockSize;

	//fast random number generator using xoshiro256+i
	//https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf	
	xso::xoshiro_4x64_plus rng;
	
	//the multiplier which will be multiplied by the 53 most significant bits of the output of the
	//random number generator to see whether a pixel should be added or not. The multiplier is
	//a uint64_t number where the least significant n bits are 1 and the remaining are 0. The pixel
	//is added when the least significant n bits (of the bit shifted 53 bits) within the rng match
	//those of the multiplier. The probability a pixel is add is then (1/2^n). Using this method,
	//knowing the amount of pixels, estimating the number of nan pixels, taking into account mindist
	//and accessible area, we can estimate a percentage chance for each pixel and set up a multiplier
	//to make that percentage happen. Doing this enables retaining only a small portion of pixel data
	//and reducing memory footprint significantly, otherwise the index of every pixel in each strata
	//would have to be stored, which would not be feasible for large rasters.
	//
	//NOTE that this method is imprecise, and has an achilles heel where if there are not very many
	//pixels of a particular strata, they may not show up at all in the final samples. For this reason,
	//the first 10000 pixels of every strata are added no matter what. If there aren't enough pixels
	//of a particular strata from the normal probability method, but there are more than 10000 pixels
	//of that strata in the raster, we may have some problems because taking from the 10000 would
	//no longer be random. Perhaps there is a way to dynamically adjust the size of the first 'thousand'
	//pixels as the raster is iterating to take this into account, without retaining too many pixels
	//to be feasible for large images.
	uint64_t multiplier = getProbabilityMultiplier(p_raster, numSamples, useMindist, access.area);

	std::vector<Index> indices;
	int64_t accessible = 0;
	int64_t notAccessible = 0;
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;

			//READ BLOCK
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			rasterBandIO(band, band.p_buffer, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
															  //read = true
			//READ ACCESS BLOCK IF NECESSARY
			if (access.used) {
				rasterBandIO(access.band, access.band.p_buffer, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false); 
			}

			//CALCULATE RAND VALUES
			for (int i = 0; i < randValIndex; i++) {
				randVals[i] = ((rng() >> 11) & multiplier) == multiplier;
			}
			randValIndex = 0;

			//PROCESS BLOCK
			switch (band.type) {
				case GDT_Int8:
					processBlock<int8_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				case GDT_UInt16:
					processBlock<uint16_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				case GDT_Int16:
					processBlock<int16_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				case GDT_UInt32:
					processBlock<uint32_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				case GDT_Int32:
					processBlock<int32_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				case GDT_Float32:
					processBlock<float>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				case GDT_Float64:
					processBlock<double>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid, accessible, notAccessible);
					break;
				default:
					throw std::runtime_error("raster pixel data type is not supported.");
			}
		}
	}

	std::shuffle(indices.begin(), indices.end(), rng);

	size_t samplesAdded = existing.used ? existing.count() : 0;
	size_t i = 0;
	while (samplesAdded < numSamples && i < indices.size()) {
		Index index = indices[i];
		i++;

		double x = GT[0] + index.x * GT[1] + index.y * GT[2];
		double y = GT[3] + index.x * GT[4] + index.y * GT[5];
		OGRPoint point = OGRPoint(x, y);

		if (mindist != 0.0 && p_layer->GetFeatureCount() != 0) {
			bool add = true;
			for (const auto &p_feature : *p_layer) {
				OGRPoint *p_point = p_feature->GetGeometryRef()->toPoint();
				if (point.Distance(p_point) < mindist) {
					add = false;
					break;
				}
			}

			if (!add) {
				continue;
			}
		}

		addPoint(&point, p_layer);
		samplesAdded++;

		if (plot) {
			xCoords.push_back(x);
			yCoords.push_back(y);
		}	
	}

	//step 10: create GDALVectorWrapper with dataset containing points
	GDALVectorWrapper *p_sampleVectorWrapper = new GDALVectorWrapper(p_samples);

	//TODO rather than first making an in-memory dataset then writing to a file afterwards,
	//just make the correct type of dataset from the get go
	if (filename != "") {
		try {
			p_sampleVectorWrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	return {{xCoords, yCoords}, p_sampleVectorWrapper, samplesAdded};
}

PYBIND11_MODULE(srs, m) {
	m.def("srs_cpp", &srs, 
		pybind11::arg("p_raster"),
		pybind11::arg("numSamples"),
		pybind11::arg("mindist"),
		pybind11::arg("p_existing").none(true),
		pybind11::arg("p_access").none(true),
		pybind11::arg("layerName"),
		pybind11::arg("buffInner"),
		pybind11::arg("buffOuter"),
		pybind11::arg("plot"),
		pybind11::arg("tempFolder"),
		pybind11::arg("filename"));

}
