/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification mapping
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

/**
 * @defgroup map
 * @ingroup stratify
 */

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "utils/raster.h"
#include "utils/helper.h"

namespace sgs {
namespace map {

/**
 * @ingroup map
 * This function maps multiple already stratified rasters into a single 
 * stratificication. Every unique combination of input band stratifications
 * corresponds to a single value in the mapped stratification.
 *
 * Based on the number of bands passed, and the number of stratifications
 * per band, multipliers are determined that will be multiplied by
 * each input band value to get the output band value.
 *
 * For example, if there were 3 stratification rasters passed, each with 5
 * possible stratum, there would be 5^3 or 125 possible mapped stratifications.
 *
 * This function can be thought of in three different sections: the setup,
 * the processing, and the finish/return. During the setup, metadata is
 * acquired for the input raster, and an output datastet is created. During
 * the procssing the input raster is iterated though, either by blocks
 * or with the entire raster in memory, the strata are determined for each
 * pixel and then written to the output dataset. During the finish/return step,
 * a GDALRasterWrapper object is created using the output dataset.
 *
 *
 * SETUP:
 * the data structures holding metadata are initialized and it is determined
 * whether the raster is a virtual raster or not, and if it is a virtual
 * raster whether it is fully in-memory or whether it must be stored on disk.
 * 
 * If the user provides an output filename, the dataset will not be a virtual dataset
 * instead it will be associated with the filename. If the user does not provide
 * an output filename then a virtual dataset driver will be used. In the case
 * of a large raster (whether or not the raster is large enough for this is calculated 
 * and passed by Python side of application), the dataset will be VRT. If the
 * package is comfortable fitting the entire raster in memory an in-memory
 * dataset will be used.
 * 
 * The input raster bands are iterated through, metadata is stored on them, and
 * bands are created for the output dataset. In the case of a VRT dataset, each
 * band is a complete dataset itself which must be added after it has been written to. 
 * In the case of a MEM dataset, the bands must be aquired from the input raster. Both MEM
 * and VRT support the AddBand function, and support bands with different types,
 * so the bands are dynamically added while iterating through the input raster
 * bands. Non virtual formats require the data types to be known at dataset initialization
 * and don't support the AddBand function, so the dataset must be created
 * after iterating through the input bands.
 *
 *
 * PROCESSING:
 * the processing section iterates through every pixel in every input band, and 
 * calculates/writes the strata to the mapped output band.
 *
 * There are two different cases -- whether the entire raster is stored in memory,
 * or whether it must be iterated though by blocks.
 *
 * If the raster is large, it is processed in blocks and splits the raster
 * into groups of blocks to be processed by multiple threads. If the raster
 * bands are in-memory, the entire raster is processed at once.For the large rasters, 
 * the processing starts out by splitting the raster into chunks depending on 
 * the number of threads. A thread is then created for each chunk. 
 * Within each thread, the blocks within it's designated chunk are iterated 
 * through and first read from the input bands, multiplied by their corresponding,
 * multipliers, then the product is written to the output band.
 *
 *
 * CLEANUP:
 * If the output dataset is a VRT dataset, the datasets which represent
 * its band (which has not yet been added as a band) must be added
 * as bands now that they are populated with data and are thus allowed
 * to be added. 
 *
 * If the dataset output band is fully in memory, it is moved to 
 * a vector from its metadata objects to be passed as a parameter
 * to the GDALRasterWrapper constructor (or not if the bands aren't in
 * memory). This GDALRasterWrapper is then returned.
 *
 *
 * @param std::vector<GDALRasterWrapper *> rasters
 * @param std::vector<std::vector<int>> bands
 * @param std::vector<std::vector<float>> strataCounts
 * @param std::string filename
 * @param bool largeRaster
 * @param in threadCount
 * @param std::string tempFolder
 * @param std::map<std::string, std::string> driverOptions,
 * @returns GDALRasterWrapper *pointer to newly created raster mapping
 */
raster::GDALRasterWrapper *map(
	std::vector<raster::GDALRasterWrapper *> rasters,
	std::vector<std::vector<int>> bands,
	std::vector<std::vector<int>> strataCounts,
	std::string filename,
	bool largeRaster,
	int threadCount,
	std::string tempFolder,
	std::map<std::string, std::string> driverOptions)
{
	GDALAllRegister();

	int height = rasters[0]->getHeight();
	int width = rasters[0]->getWidth();
	double *geotransform = rasters[0]->getGeotransform();
	std::string projection = (rasters[0]->getDataset()->GetProjectionRef());
	if (projection == "") {
		throw std::runtime_error("could not get projection from the first raster argument.");
	}

	OGRSpatialReference srs;
	srs.importFromWkt(projection.c_str());

	for (size_t i = 1; i < rasters.size(); i++) {
		if (rasters[i]->getHeight() != height) {
			std::string err = "raster with index " + std::to_string(i) + " has a different height from the raster at index 0.";
			throw std::runtime_error(err);
		}

		if (rasters[i]->getWidth() != width) {
			std::string err = "raster with index " + std::to_string(i) + " has a different width from the raster at index 0.";
			throw std::runtime_error(err);
		}

		double *checkGeotransform = rasters[i]->getGeotransform();
		for (int j = 0; j < 6; j++) {
			if (geotransform[i] != checkGeotransform[i]) {
				std::string err = "raster with index " + std::to_string(i) + " has a different geotransform from the raster at index 0.";
				throw std::runtime_error(err);
			}
		}

		OGRSpatialReference checkSrs;
		checkSrs.importFromWkt(rasters[i]->getDataset()->GetProjectionRef());
		if (!srs.IsSame(&checkSrs)) {
			std::string err = "raster with index " + std::to_string(i) + " has a different projection from the raster at index 0.";
			throw std::runtime_error(err);
		}
	}
	
	std::vector<helper::RasterBandMetaData> stratBands;
	std::vector<int> numStrataPerBand;	
	helper::RasterBandMetaData mapBand;
	std::vector<helper::VRTBandDatasetInfo> VRTBandInfo(1);

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";

	std::vector<std::mutex> stratDatasetMutexes(rasters.size());
	std::mutex mapBandMutex;

	std::string driver;
	GDALDataset *p_dataset = nullptr;
	if (isMEMDataset || isVRTDataset) {
		std::string driver = isMEMDataset ? "MEM" : "VRT";
		p_dataset = helper::createVirtualDataset(driver, width, height, geotransform, projection);
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();

		if (extension == ".tif") {
			driver = "Gtiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now");
		}
	}

	//step 1 iterate through bands populating rasterBands and bandStratMultiplier objects
	std::vector<size_t> multipliers(1, 1);
	for (size_t i = 0; i < rasters.size(); i++) {
		raster::GDALRasterWrapper *p_raster = rasters[i];

		for (size_t j = 0; j < bands[i].size(); j++) {
			int band = bands[i][j];
			int strataCount = strataCounts[i][j];
			numStrataPerBand.push_back(strataCount);

			stratBand = helper::rasterBandMetaData;

			GDALRasterBand *p_band = p_raster->getRasterBand(band);
			stratBand.p_band = p_band;
			stratBand.type = p_raster->getRasterBandType(band);
			stratBand.size = p_raster->getRasterBandTypeSize(band);
			stratBand.p_buffer = largeRaster ? nullptr : p_raster->getRasterBandBuffer(band);
			stratBand.nan = p_band->GetNoDataValue();
			stratBand.p_mutex = &stratDatasetMutexes[i];
			p_band->GetBlockSize(&stratBand.xBlockSize, &stratBand.yBlockSize);
			stratBands.push_back(stratBand);

			helper::printTypeWarningsForInt32Conversion(p_stratBand->type);
			multipliers.push_back(multipliers.back() * strataCount);
		}
	}

	size_t bandCount = stratBands.size();
	size_t maxStrata = multipliers.back();
	multipliers.pop_back();
	helper::setStratBandTypeAndSize(maxStrata, &mapBand.type, &mapBand.size);
	mapBand.name = "strat_map";
	mapBand.xBlockSize = stratBands[0].xBlockSize;
	mapBand.yBlockSize = stratBands[0].yBlockSize;
	mapBand.p_mutex = &mapBandMutex;

	if (isMEMDataset) {
		helper::addBandToMEMDataset(p_dataset, mapBand);
	}
	else if (isVRTDataset) {
		helper::createVRTBandDataset(
			p_dataset,
			mapBand,
			tempFolder,
			"map",
			VRTBandInfo,
			driverOptions
		);
	}
	else {
		bool useTiles = mapBand.xBlockSize != width &&
				mapBand.yBlockSize != height;

		mapBand.p_buffer = !largeRaster ?
			VSIMalloc3(height, width, mapBand.size) :
			nullptr;

		p_dataset = helper::createDataset(
			filename,
			driver,
			width,
			height,
			geotransform,
			projection,
			&mapBand,
			1,
			useTiles,
			driverOptions
		);
	}

	if (largeRaster) {
		pybind11::gil_scoped_acquire acquire;
		boost::asio::thread_pool pool(threadCount);

		int xBlockSize = stratBands[0].xBlockSize;
		int yBlockSize = stratBands[0].yBlockSize;

		int xBlocks = (width + xBlockSize - 1) / xBlockSize;
		int yBlocks = (height + yBlockSize - 1) / yBlockSize;
		int chunkSize = yBlocks / threadCount;

		for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
			int yBlockEnd = std::min(yBlockStart + chunkSize, yBlocks);

			boost::asio::post(pool, [
				bandCount,
				xBlockSize,
				yBlockSize,
				yBlockStart,
				yBlockEnd,
				xBlocks,
				&stratBands,
				&mapBand,
				&multipliers	
			] {
				std::vector<void *> stratBuffers(bandCount);
				for (size_t band = 0; band < bandCount; band++) {
					stratBuffers[band] = VSIMalloc3(xBlockSize, yBlockSize, stratBands[band].size);
				}
				void *p_mapBuffer = VSIMalloc3(xBlockSize, yBlockSize, mapBand.size);

				//cast to int so there is less data type conversion while processign
				std::vector<int> intNoDataValues(bandCount);
				for (size_t band = 0; band < bandCount; band++) {
					intNoDataValues[band] = static_cast<int>(stratBands[band].nan);
				}

				for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
					for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
						int xValid, yValid;
						stratBands[0].p_mutex->lock();
						stratBands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
						stratBands[0].p_mutex->unlock();

						for (size_t band = 0; band < bandCount; band++) {
							rasterBandIO(
								stratBands[band],
								stratBuffers[band],
								xBlockSize,
								yBlockSize,
								xBlock,
								yBlock,
								xValid,
								yValid,
								true //read = true
							);
						}

						for (int y = 0; y < yValid;  y++) {
							for (int x = 0; x < xValid; x++) {
								size_t index = static_cast<size_t>(x + y * xBlockSize);
								bool isNan = false;
								int mappedStrat = 0;

								for (size_t band = 0; band < bandCount; band++) {
									int strat = helper::getPixelValueDependingOnType<int>(stratBands[band].type, stratBuffers[band], index);
									isNan = strat == intNoDataValues[band];
	
									if (!isNan) {
										if (strat >= numStrataPerBand[band]) {
											std::string errmsg = "the num_strata indicated for band " + dataBands[band].p_band->GetDescription() " is less than or equal to one of the values in that band.";
											throw std::runtime_error(errmsg);
										}

										if (strat < 0) {
											std::string errmsg = "a negative strata value of " + std::to_string(strat) + " was found in band " + dataBands[band].p_band->GetDescription() ", and is not marked as a nodata value.";
											throw std::runtime_error(errmsg);
										}

										mappedStrat += strat * multipliers[band];
									}
									else {
										break;
									}
								}

								helper::setStrataPixelDependingOnType(mapBand.type, p_mapBuffer, index, isNan, mappedStrat);
							}
						}

						helper::rasterBandIO(
							mapBand,
							p_mapBuffer,
							xBlockSize,
							yBlockSize,
							xBlock,
							yBlock,
							xValid,
							yValid,
							false //read = false
						);
					}
				}

				for (size_t band = 0; band < bandCount; band++) {
					VSIFree(stratBuffers[band]);
				}
				VSIFree(p_mapBuffer);
			});
		}
		
		pool.join();
		pybind11::gil_scoped_release release;
	}
	else {
		std::vector<int> intNoDataValues(bandCount);
		for (size_t band = 0; band < bandCount; band++) {
			intNoDataValues[band] = static_cast<int>(stratBands[band].nan);
		}

		size_t pixelCount = static_cast<size_t>(height) * static_cast<size_t>(width);
		for (size_t index = 0; index < pixelCount;  index++) {
			bool isNan = false;
			int mappedStrat = 0;

			for (size_t band = 0; (band < bandCount) && !isNan; band++) {
				int strat = helper::getPixelValueDependingOnType<int>(stratBands[band].type, stratBands[band].p_buffer, index);
				isNan = strat == intNoDataValues[band];
				if (!isNan) {
					mappedStrat += strat * multipliers[band];
				}
			}

			helper::setStrataPixelDependingOnType(mapBand.type, mapBand.p_buffer, index, isNan, mappedStrat);
		}

		if (!isVRTDataset && !isMEMDataset) {
			CPLErr err = mapBand.p_band->RasterIO(
				GF_Write,
				0,
				0,
				width,
				height,
				mapBand.p_buffer,
				width,
				height,
				mapBand.type,
				0,
				0
			);
			if (err) {
				throw std::runtime_error("error writing band to file.");
			}
		}
	}

	//close and add VRT sub dataset as band
	if (isVRTDataset) {
		GDALClose(VRTBandInfo[0].p_dataset);
		helper::addBandToVRTDataset(p_dataset, mapBand, VRTBandInfo[0]);
	}

	return largeRaster ?
		new raster::GDALRasterWrapper(p_dataset) :
		new raster::GDALRasterWrapper(p_dataset, {mapBand.p_buffer});
}

} //namespace map
} //namespace sgs
