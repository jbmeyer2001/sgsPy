/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of random sampling
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

/**
 * @defgroup srs srs
 * @ingroup sample
 */

#include <iostream>
#include <random>
#include <unordered_map>

#include <xoshiro.h>

#include "utils/access.h"
#include "utils/existing.h"
#include "utils/helper.h"
#include "utils/raster.h"
#include "utils/vector.h"

namespace sgs {
namespace srs {

/**
 * @ingroup srs
 *
 * This function generates random index values to sample. This method is fast
 * in many circumstances, because it does not require the entire raster to be read.
 * 
 * However, in cases where there are a very large number of samples, or a low percentage
 * of the raster is sampleable, it may be slower than reading through the entire raster.
 * This is because calling the RasterIO function incurs significant overhead, especially
 * in random access patterns.
 *
 * @param helper::RasterBandMetaData& band
 * @param int width
 * @param int height
 * @param int numSamples
 * @param access::Access& access
 * @param existing::Existing& existing
 * @param std::unordered_set<helper::Index>& index
 * @param xso::xoshiro_4x64_plus rng
 * @returns bool
 */
template <typename T>
inline bool
getRandomIndices(
	helper::RasterBandMetaData& band,
	int width,
	int height,
	int numSamples,
	access::Access& access,
	existing::Existing& existing,
	std::vector<helper::Index>& indices,
	xso::xoshiro_4x64_plus rng)
{
	T nan = static_cast<T>(band.nan);

	uint64_t width64 = static_cast<int64_t>(width);
	uint64_t height64 = static_cast<int64_t>(height);

	//get mask for rng output using bit twiddling
	uint64_t maxIndex = width64 * height64 - 1;
	uint64_t mask = maxIndex;
	mask |= mask >> 1;
	mask |= mask >> 2;
	mask |= mask >> 4;
	mask |= mask >> 8;
	mask |= mask >> 16;
	mask |= mask >> 32;

	int iterations = 0;
	int maxIterations = numSamples * 11;
	std::unordered_set<uint64_t> indexSet;
	while (iterations < maxIterations && indices.size() < static_cast<size_t>(numSamples)) {
		//generate a random valid index
		//
		//bit shift by 11 because lower 11 bits are not random enough in this (fast) xoshiro configuration
		uint64_t index = (rng() >> 11) & mask; 
		while (index > maxIndex) {
			index = (rng() >> 11) & mask;
		}
		int x = static_cast<int>(index % width64);
		int y = static_cast<int>(index / width64);
		
		//read the value from that index
		T val;
		CPLErr err = band.p_band->RasterIO(GF_Read, x, y, 1, 1, &val, 1, 1, band.type, 0, 0);
		if (err) {
			throw std::runtime_error("error reading pixel from raster band using GDALRasterBand::RasterIO().");
		}

		//check if val is nan
		bool isNan = std::isnan(val) || val == nan;
		if (isNan) {
			iterations++;
			continue;
		}

		//check if pixel is accessible
		bool accessible = !access.used;
		if (access.used) {
			int8_t accessVal;
			CPLErr err = access.band.p_band->RasterIO(GF_Read, x, y, 1, 1, &accessVal, 1, 1, access.band.type, 0, 0);
			if (err) {
				throw std::runtime_error("error reading pixel from access raster band using GDALRasterBand::RasterIO().");
			}
			accessible = accessVal != 1;
		}
		if (!accessible) {
			iterations++;
			continue;
		}

		//check if pixel is already sampled
		bool alreadySampled = existing.used && existing.containsIndex(x, y);
		alreadySampled |= indexSet.find(index) != indexSet.end();
		if (alreadySampled) {
			iterations++;
			continue;
		}

		//add index to indices map if the pixel is not nan, accessible, and not already sampled
		indexSet.insert(index);
		indices.push_back(helper::Index(x, y));
		iterations++;	
	}

	return indices.size() == static_cast<size_t>(numSamples);
}

/**
 * @ingroup srs
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
	helper::RasterBandMetaData& band,
	access::Access& access,
	existing::Existing& existing,
	std::vector<helper::Index>& indices,
	helper::RandValController& rand,
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
			//get val
			T val = helper::getPixelValueDependingOnType<T>(band.type, band.p_buffer, blockIndex);

			//check nan
			bool isNan = val == nan || std::isnan(val);
			if (isNan) {
				blockIndex++;
				continue;
			}

			//check access
			if (access.used && p_access[blockIndex] == 1) {
				blockIndex++;
				continue;
			}

			helper::Index index = {x + xBlock * band.xBlockSize, y + yBlock * band.yBlockSize};
			
			//check existing
			if (existing.used && existing.containsIndex(index.x, index.y)) {
				blockIndex++;
				continue;
			}

			//check rand val
			if (!rand.next()) {
				blockIndex++;
				continue;
			}

			//add index to indices
			indices.push_back(index);

			//increment within-block index
			blockIndex++;
		}
	}
}

/**
 * @ingroup srs
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
 * Then, the raster is processed in one of two ways. One using a random access
 * strategy where random indexes are calculated and checked for validity, and 
 * another where the entire raster is read. The decision between these two
 * is calculated to minimize the number of blocks read into memory, as this
 * is the main bottleneck in processing time.
 *
 * Once all possible pixels have been selected, there may be extra
 * indices in the indicies vector. Because simply sampling the first
 * few we need would might NOT be in a random order, the indices are first shuffled.
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
 * @returns std::tuple<std::vector<std::vector<double>>, GDALVectorWrapper *, size_t>
 */
std::tuple<std::vector<std::vector<double>>, vector::GDALVectorWrapper *, size_t> 
srs(
	raster::GDALRasterWrapper *p_raster,
	size_t numSamples,
	double mindist,
	vector::GDALVectorWrapper *p_existing,
	vector::GDALVectorWrapper *p_access,
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
	helper::RasterBandMetaData band;
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

	vector::GDALVectorWrapper *p_wrapper = new vector::GDALVectorWrapper(p_samples, std::string(p_raster->getDataset()->GetProjectionRef()));
	OGRLayer *p_layer = p_samples->CreateLayer("samples", p_wrapper->getSRS(), wkbPoint, nullptr);
	if (!p_layer) {
		throw std::runtime_error("unable to create output dataset layer.");
	}

	//generate access structure
	access::Access access(
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
	existing::Existing existing(
		p_existing,
		p_raster,
		GT,
		p_raster->getWidth(),
		p_layer,
		plot,
		xCoords,
		yCoords
	);

	int xBlocks = (width + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (height + band.yBlockSize - 1) / band.yBlockSize;
	std::vector<helper::Index> indices;

	std::vector<bool> randVals(band.xBlockSize * band.yBlockSize);
	int randValIndex = band.xBlockSize * band.yBlockSize;

	//fast random number generator using xoshiro256++
	//https://vigna.di.unimi.it/ftp/papers/ScrambledLinear.pdf	
	xso::xoshiro_4x64_plus rng;

	// when reading pixels or blocks into memory using GDAL, the whole block is always read into memory.
	// This reading of blocks is a large portion of the runtime for the processing of raster images.
	//
	// When using the Simple Random Sampling (srs) method, there is no requirement that the pixels
	// have any relation to each other or represent any kind of distribution in the overall raster.
	// This means that there is no requirement to read through the entire raster image.
	//
	// In some cases, it may be faster to select pixels at random from the image, hoping they are
	//  1. not nan
	//  2. acessible
	//  3. not already sampled
	//
	// However, there are cases when more blocks will be read into memory by randomly selecting 
	// pixels compared to reading through the whole raster. This is because random access patterns
	// have a VERY low hit rate in the cache, leading to significant thrashing -- we have to assume
	// every time a pixel is read a block is read into the cache whereas when the entire raster is
	// read sequentially, a new block must be read into the cache ONLY when the first pixel in that
	// block is looked at.
	//
	// In the following, we attemt to estimate which method will result in the smaller number of 
	// blocks read into memory and use that to generate random pixels, assuming that half of the
	// raster is nan.
	//
	// worth noting -- if not enough pixels are able to be determined using the random access method,
	// the full raster read method is then used.

	//total blocks is yBlocks * xBlocks	
	size_t numBlocks = static_cast<size_t>(xBlocks) * static_cast<size_t>(yBlocks);

	//desired samples is x3 if using mindist
	size_t desiredSamples = static_cast<size_t>(useMindist ? numSamples * 3 : numSamples);

	//total number of pixels in the raster
	size_t allPixels = static_cast<size_t>(width) * static_cast<size_t>(height);

	//total area is width * height * pixelWidth * pixelHeight == allPixels * pixelWidth * pixelHeight
	double totalArea = static_cast<double>(p_raster->getPixelHeight()) * 
		           static_cast<double>(p_raster->getPixelWidth()) * 
			   static_cast<double>(allPixels);

	//accessible pixels is the proportion of pixels which are within the accessible area times the total pixel count
	double accessiblePixels = allPixels * (access.used ? (access.area / totalArea) : 1.0);

	//valid pixels divides the accessible pixels by 2 (assuming half nan), then subtracts existing samples
	size_t validPixels = static_cast<size_t>(accessiblePixels / 2.0) - (existing.used ? existing.samples.size() : 0);
	
	size_t randomAccessBlocksRequired = 0;

        while (desiredSamples > 0 && randomAccessBlocksRequired < numBlocks) {
		//the likelihood of guessing a valid pixel is equal to (remainingValidPixles / allPixels)
		//
		//therefore it is estimated there will be a total number of 1 / (remainingValidPixels / allPixels)
		//guesses or random block accesses.
		//
		//This is equivalent to allPixels / validPixels.
		//
		//Then, since we sampled 1 pixel we reduce the number of valid remaining pixels by 1.

		randomAccessBlocksRequired += (allPixels / validPixels);
		validPixels--;
		desiredSamples--;
	}

	//set max number of random accesses as the estimated required number * 1.25
	size_t maxRandomAccessBlocks = randomAccessBlocksRequired + randomAccessBlocksRequired / 4;

	//if the max estimated number of blocks read into memory under a random access strategy is less than
	//reading the whole raster, try the random access strategy capping the number of accesses at this max value.
	//
	//Then, only read the entire raster if not enough pixels were read by the random strategy
	bool haveEnoughSamples = false;
	if (maxRandomAccessBlocks < numBlocks) {
		desiredSamples = static_cast<size_t>(useMindist ? numSamples * 3 : numSamples);

		switch (band.type) {
			case GDT_Int8:
				haveEnoughSamples = getRandomIndices<int8_t>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			case GDT_UInt16:
				haveEnoughSamples = getRandomIndices<uint16_t>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			case GDT_Int16:
				haveEnoughSamples = getRandomIndices<int16_t>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			case GDT_UInt32:
			 	haveEnoughSamples = getRandomIndices<uint32_t>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			case GDT_Int32:
				haveEnoughSamples = getRandomIndices<int32_t>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			case GDT_Float32:
				haveEnoughSamples = getRandomIndices<float>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			case GDT_Float64:
				haveEnoughSamples = getRandomIndices<double>(band, width, height, desiredSamples, access, existing, indices, rng);
				break;
			default:
				throw std::runtime_error("raster pixel data type not supported.");
		}
	}

	if (!haveEnoughSamples) {
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
		uint64_t multiplier = helper::getProbabilityMultiplier(
			width, 
			height, 
			p_raster->getPixelWidth(), 
			p_raster->getPixelHeight(), 
			8, 
			numSamples, 
			useMindist, 
			access.area
		);
	
		helper::RandValController rand(band.xBlockSize, band.yBlockSize, multiplier, &rng);

		for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
			for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
				int xValid, yValid;
	
				//read block
				band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
				helper::rasterBandIO(band, band.p_buffer, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false);
																//read = true
				//read access block if necessary
				if (access.used) {
					helper::rasterBandIO(access.band, access.band.p_buffer, band.xBlockSize, band.yBlockSize, xBlock, yBlock, xValid, yValid, true, false); 
				}
	
				//calculate random values
				rand.calculateRandValues();
	
				//process block
				switch (band.type) {
					case GDT_Int8:
						processBlock<int8_t>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					case GDT_UInt16:
						processBlock<uint16_t>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					case GDT_Int16:
						processBlock<int16_t>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					case GDT_UInt32:
						processBlock<uint32_t>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					case GDT_Int32:
						processBlock<int32_t>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					case GDT_Float32:
						processBlock<float>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					case GDT_Float64:
						processBlock<double>(band, access, existing, indices, rand, xBlock, yBlock, xValid, yValid);
						break;
					default:
						throw std::runtime_error("raster pixel data type is not supported.");
				}
			}
		}	
	}

	std::shuffle(indices.begin(), indices.end(), rng);
	
	size_t samplesAdded = existing.used ? existing.count() : 0;
	size_t i = 0;
	helper::NeighborMap neighbor_map;
	double mindist_sq = mindist * mindist;
	
	helper::Field fieldExistingFalse("existing", 0);
	while (samplesAdded < numSamples && i < indices.size()) {
		helper::Index index = indices[i];
		bool valid = true;
		const auto [x, y] = helper::sample_to_point(GT, index);
	
		if (useMindist) {
			valid = helper::is_valid_sample(x, y, neighbor_map, mindist, mindist_sq);
		}
	
		if (valid) {
			OGRPoint point = OGRPoint(x, y);
	
			existing.used ?
				helper::addPoint(&point, p_layer, &fieldExistingFalse) :
				helper::addPoint(&point, p_layer);
	
			samplesAdded++;
	
			if (plot) {
				xCoords.push_back(x);
				yCoords.push_back(y);
			}
		}
		i++;
	}


	if (filename != "") {
		try {
			p_wrapper->write(filename);
		}
		catch (const std::exception& e) {
			std::cout << "Exception thrown trying to write file: " << e.what() << std::endl;
		}
	}

	//free allocated memory
	VSIFree(band.p_buffer);
	if (access.used) {
		VSIFree(access.band.p_buffer);
	}

	return {{xCoords, yCoords}, p_wrapper, samplesAdded};
}

} //namespace srs
} //namespace sgs
