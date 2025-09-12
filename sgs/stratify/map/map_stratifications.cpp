/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification mapping
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "raster.h"
#include "helper.h"

/**
 * This function maps multiple stratifications.
 *
 * Based on the number of bands passed, and the number of stratifications
 * per band, multipliers are determined such that every unique combination
 * of stratifications within the passed stratification rasters
 * corrosponds to a single stratification in the output mapping.
 *
 * For example, if there were 3 stratification rasters passed, each with 5
 * possible stratum, there would be 5^3 or 125 possible mapped stratifications.
 *
 * A new GDALRasterWrappper object is created and initialized using the
 * mapped raster, and the mapped raster is written to disk if a filename is given.
 *
 * NOTE: stratifications have the data type 'float' in order to accomodate nan values
 *
 * @param std::vector<GDALRasterWrapper *> rasters list
 * @param std::vector<std::vector<int>> bands list
 * @param std::vector<std::vector<float>> number of stratums list
 * @param std::string filename the filename to write to or "" if file shouldn't be written
 * @returns GDALRasterWrapper *pointer to newly created raster mapping
 */
GDALRasterWrapper *mapStratifications(
	std::vector<GDALRasterWrapper *> rasters,
	std::vector<std::vector<int>> bands,
	std::vector<std::vector<int>> strataCounts,
	std::string filename,
	bool largeRaster,
	int threads,
	std::string tempFolder,
	std::map<std::string, std::string> driverOptions)
{
	GDALAllRegister();

	//define useful variables
	int height = rasters[0]->getHeight();
	int width = rasters[0]->getWidth();
	double *geotransform = rasters[0]->getGeotransform();
	std::string projection = std::string(rasters[0]->getDataset()->GetProjectionRef());

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
	}
	
	std::vector<RasterBandMetaData> stratBands;
	RasterBandMetaData mapBand;
	std::vector<VRTBandDatasetInfo> VRTBandInfo(1);

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";

	std::vector<std::mutex> stratDatasetMutexes(rasters.size());
	std::mutex mapBandMutex;

	std::string driver;
	GDALDataset *p_dataset = nullptr;
	if (isMEMDataset || isVRTDataset) {
		std::string driver = isMEMDataset ? "MEM" : "VRT";
		p_dataset = createVirtualDataset(driver, width, height, geotransform, projection);
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
		GDALRasterWrapper *p_raster = rasters[i];

		for (size_t j = 0; j < bands[i].size(); j++) {
			stratBands.resize(stratBands.size() + 1);
			RasterBandMetaData *p_stratBand = &stratBands.back();
			
			int band = bands[i][j];
			int strataCount = strataCounts[i][j];

			GDALRasterBand *p_band = p_raster->getRasterBand(band);
			p_stratBand->p_band = p_band;
			p_stratBand->type = p_raster->getRasterBandType(band);
			p_stratBand->size = p_raster->getRasterBandTypeSize(band);
			p_stratBand->p_buffer = largeRaster ? nullptr : p_raster->getRasterBandBuffer(band);
			p_stratBand->nan = p_band->GetNoDataValue();
			p_stratBand->p_mutex = &stratDatasetMutexes[i];
			p_band->GetBlockSize(&p_stratBand->xBlockSize, &p_stratBand->yBlockSize);

			printTypeWarningsForInt32Conversion(p_stratBand->type);
			multipliers.push_back(multipliers.back() * strataCount);
		}
	}

	size_t bandCount = stratBands.size();
	size_t maxStrata = multipliers.back();
	multipliers.pop_back();
	setStratBandTypeAndSize(maxStrata, &mapBand.type, &mapBand.size);
	mapBand.name = "strat_map";
	mapBand.xBlockSize = stratBands[0].xBlockSize;
	mapBand.yBlockSize = stratBands[0].yBlockSize;
	mapBand.p_mutex = &mapBandMutex;

	if (isMEMDataset) {
		addBandToMEMDataset(p_dataset, mapBand);
	}
	else if (isVRTDataset) {
		createVRTBandDataset(
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

		p_dataset = createDataset(
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
		boost::asio::thread_pool pool(1); //TODO switch back
		//boost::asio::thread_pool pool(threads);

		int xBlockSize = stratBands[0].xBlockSize;
		int yBlockSize = stratBands[0].yBlockSize;

		int xBlocks = (width + xBlockSize - 1) / xBlockSize;
		int yBlocks = (height + yBlockSize - 1) / yBlockSize;
		int chunkSize = yBlocks / threads;

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
									int strat = getPixelValueDependingOnType<int>(stratBands[band].type, stratBuffers[band], index);
									isNan = strat == intNoDataValues[band];

									if (!isNan) {
										mappedStrat += strat * multipliers[band];
									}
									else {
										break;
									}
								}

								setStrataPixelDependingOnType(mapBand.type, p_mapBuffer, index, isNan, mappedStrat);
							}
						}

						rasterBandIO(
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
		for (size_t band; band < bandCount; band++) {
			intNoDataValues[band] = static_cast<int>(stratBands[band].nan);
		}

		size_t pixelCount = static_cast<size_t>(height) * static_cast<size_t>(width);
		for (size_t index = 0; index < pixelCount;  index++) {
			bool isNan = false;
			int mappedStrat = 0;

			for (size_t band = 0; (band < bandCount) && !isNan; band++) {
				int strat = getPixelValueDependingOnType<int>(stratBands[band].type, stratBands[band].p_buffer, index);
				isNan = strat == intNoDataValues[band];
				if (!isNan) {
					mappedStrat += strat * multipliers[band];
				}
			}

			setStrataPixelDependingOnType(mapBand.type, mapBand.p_buffer, index, isNan, mappedStrat);
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
		addBandToVRTDataset(p_dataset, mapBand, VRTBandInfo[0]);
	}

	return largeRaster ?
		new GDALRasterWrapper(p_dataset) :
		new GDALRasterWrapper(p_dataset, {mapBand.p_buffer});
}

PYBIND11_MODULE(map_stratifications, m) {
	m.def("map_stratifications_cpp", &mapStratifications);
}
