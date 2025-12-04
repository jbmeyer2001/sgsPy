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

#include "utils/access.h"
#include "utils/existing.h"
#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

namespace srs {

/**
 * This is a helper function for processing a block of the raster. For each
 * pixel in the block: 
 * The value is checked, and not added if it is a nanvalue. 
 * The pixel is checked to ensure it is within an accessible area.
 * The pixel is checked to ensure it hasn't already been added as a pre-existing sample point.
 * The next rng value is checked to see whether it is one of the chosen pixels to be added.
 *
 * @param RasterBandMetaData& band
 * @param Access& access
 * @param Existing& existing
 * @param std::vector<Index>& indices
 * @param std::vector<Index>& randVals,
 * @param int& randValIndex
 * @param int xBlock
 * @param int yBlock
 * @param int xValid
 * @param int yValid 
 */
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
	int yValid)
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
			if (access.used && p_access[blockIndex] == 1) {
				blockIndex++;
				continue;
			}

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
 * This function uses random sampling to determine the location
 * of sample plots given a raster image.
 *
 * First, metadata is acquired on the first raster band, which is 
 * to be read to check and ensure samples don't occur over nodata pixels.
 *
 * Next, and output vector dataset is created as an in-memory dataset.
 * If the user specifies a filename, this in-memory dataset will be 
 * written to disk in a different format after all points have been added.
 *
 * An Access struct is created, which creates a raster dataset containing
 * a rasterized version of access buffers. This raster will be 1 over
 * accessible areas. In the case where there is no access vector given,
 * the structs 'used' member will be false and no processing or rasterization
 * will be done.
 *
 * An Existing struct is created, which retains information on already existing
 * sample points passed in the form of a vector dataset. The points are iterated
 * through and added to the output dataset. The points are also added to a set,
 * and during iteration the indexes of every pixel will be checked against this set
 * to ensure there are no duplicate pixels. In the case whre there is no existing
 * vector given, the structs 'used' member will be false and no processing
 * will be done.
 *
 * Next, a rng() function is created usign the xoshiro library, the specific
 * randm number generator is the xoshrio256++
 * https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf	
 *
 * The impetus behind usign the rng() function to determine which pixels
 * should be added DURING iteration, rather than afterwards, is it removes the
 * necessity of storing every available pixel, which quickly becomes extrordinarily
 * inefficient for large rasters. Rather, for pixels which are accessible, not nan,
 * and not already existing, there is a pre-determined percentage chance to be stored
 * which uses this random number generator. An over-estimation for the percentage
 * chance is made, because it is better to have too many than not enough possible options
 * to sample from. This over-estimation might result in the storage of 2x-3x extra pixels
 * rather than the many orders of magnitude extra storage of adding all pixels. The
 * calculation for this percentage is done and explained in detail in the
 * getProbabilityMultiplier() function.
 *
 * The raster is then processed in blocks, each block is read into memory,
 * and potentially the block of the access raster is read into memory as well.
 * The processBlock() funciton is called depending on the data type of the 
 * raster, to add the available / chosen pixel indices.
 *
 * Once the entire raster has been iterated through, there may be extra
 * indices in the indicies vector. Because simply sampling the first
 * few we need would NOT be random, the indices are first shuffled.
 * After being shuffled, the indexes are added to the output dataset
 * as samples if they don't occur within mindist if an already existing pixel.
 *
 * @param GDALRasterWrapper *p_raster
 * @param size_t numSamples
 * @param double mindist
 * @param GDALVectorWrapper *p_existing
 * @param GDALVectorWrapper *p_access
 * @param std::string layerName
 * @param double buffInner
 * @param double buffOuter
 * @param bool plot
 * @param std::string tempFolder
 * @param std::string filename
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

	//fast random number generator using xoshiro256++
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
	//and reducing memory footprint significantly, otherwise the index of every pixel
	//would have to be stored, which would not be feasible for large rasters.
	uint64_t multiplier = getProbabilityMultiplier(
		width, 
		height, 
		p_raster->getPixelWidth(), 
		p_raster->getPixelHeight(), 
		4, 
		numSamples, 
		useMindist, 
		access.area
	);

	std::vector<Index> indices;
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
					processBlock<int8_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
					break;
				case GDT_UInt16:
					processBlock<uint16_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
					break;
				case GDT_Int16:
					processBlock<int16_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
					break;
				case GDT_UInt32:
					processBlock<uint32_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
					break;
				case GDT_Int32:
					processBlock<int32_t>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
					break;
				case GDT_Float32:
					processBlock<float>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
					break;
				case GDT_Float64:
					processBlock<double>(band, access, existing, indices, randVals, randValIndex, xBlock, yBlock, xValid, yValid);
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

} //namespace srs
