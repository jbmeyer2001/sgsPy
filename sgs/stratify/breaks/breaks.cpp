/******************************************************************************
 *
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using breaks
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "raster.h"
#include "helper.h"

/**
 * This is a helper function for processing a pixel of data
 * when a mapped stratification is being created.
 *
 * First, the value is read in as a double, and it is determined
 * whether the pixel is a nan pixel or not. The mapNan boolean 
 * is updated in addition to the isNan boolean, to ensure that
 * if one band within the raster is nan at a certain pixel
 * then the mapped raster (but not necessarily all output rasters)
 * is also nan at that pixel.
 *
 * Then, if it isn't a nan pixel the lower bound of the value 
 * within the vector of break values is found. For example,
 * if the value was 3 and the breaks vector was [2, 4, 6], the
 * lower bound would be 1, which is the index of 4, the first
 * value larger than 3 in the breaks vector. This lower bound
 * is the strata. This strata (or the nan value) is then written 
 * with the appropriate type to the strat raster band.
 */
inline void processMapPixel(
	size_t index,
	RasterBandMetaData& dataBand,
	void * p_dataBuffer,
	RasterBandMetaData& stratBand,
	void * p_stratBuffer,
	std::vector<double>& bandBreaks,
	size_t multiplier,
	bool& mapNan,
	size_t& mapStrat)
{
	double val = getPixelValueDependingOnType<double>(dataBand.type, p_dataBuffer, index);
	bool isNan = std::isnan(val) || (double)val == dataBand.nan;
	mapNan |= isNan;
		
	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(bandBreaks.begin(), bandBreaks.end(), val);
		strat = (it == bandBreaks.end()) ? bandBreaks.size() : std::distance(bandBreaks.begin(), it);
	}

	setStrataPixelDependingOnType(stratBand.type, p_stratBuffer, index, isNan, strat);

	if (!mapNan) {
		mapStrat += strat * multiplier;
	}
}	

/**
 * This is a helper function for processing a pixel of data.
 *
 * First, the value is read in as a double, and it is determined
 * whether the pixel is a nan pixel or not.
 *
 * Then, if it isn't a nan pixel the lower bound of the value 
 * within the vector of break values is found. For example,
 * if the value was 3 and the breaks vector was [2, 4, 6], the
 * lower bound would be 1, which is the index of 4, the first
 * value larger than 3 in the breaks vector. This lower bound
 * is the strata. This strata (or the nan value) is then written 
 * with the appropriate type to the strat raster band.
 */
inline void
processPixel(
	size_t index,
	void *p_data,
	RasterBandMetaData *p_dataBand,
	void *p_strat,
	RasterBandMetaData *p_stratBand,
	std::vector<double>& bandBreaks)
{
	double val = getPixelValueDependingOnType<double>(p_dataBand->type, p_data, index);
	bool isNan = std::isnan(val) || val == p_dataBand->nan;

	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(bandBreaks.begin(), bandBreaks.end(), val);
		strat = (it == bandBreaks.end()) ? bandBreaks.size() : std::distance(bandBreaks.begin(), it);
	}

	setStrataPixelDependingOnType(p_stratBand->type, p_strat, index, isNan, strat);
}

/**
 * This function stratifies a given raster using user-defined breaks.
 * The breaks are provided as a vector of doubles for each band specified
 * in the input dataset.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map variable to combine the stratification of
 * multiple raster bands.
 *
 * The function can be thought of in three different sections: the setup,
 * the processing, and the finish/return. During the setup, metadata is aquired
 * for the input raster, and an output dataset is created which depends
 * on user-given parameters and the input raster. During the processing
 * the input raster is iterated through, either by blocks or with the
 * entire raster in memory, the strata are determined for each pixel and
 * then written to the output dataset. During the finish/return step, a
 * GDALRasterWrapper object is created using the output dataset. 
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
 * calculates/writes the strata to the corresponding output band.
 *
 * There are four different cases dealing with whether or not the entire raster
 * band is allocated in memory (the largeRaster variable is false), and whether
 * or not the values of each band should be mapped to an extra output raster
 * band.
 *
 * If the raster is large, it is processed in blocks and splits the raster
 * into groups of blocks to be processed by multiple threads. If the raster
 * bands are in-memory, the entire raster is processed at once. The mapped
 * rasters store information on an extra output raster band, the output
 * values of which are determined as a function of all other output raster
 * bands. The multipliers vector stores the information for this.
 *
 * For the large rasters, the processing starts out by splitting the raster
 * into chunks depending on the number of threads. A thread is then created
 * for each chunk. Within each thread, the blocks within it's designated
 * chunk are iterated through and first read from the input bands, processed,
 * then written to the output bands. In the case of a mapped raster all of the
 * bands are iterated alongside eachother so that the intermediate mapping
 * calculations don't have to be written then read again. In the case of a non
 * mapped raster, each band is processed sequentially.
 *
 * CLEANUP:
 * If the output dataset is a VRT dataset, the datasets which represent
 * its bands (that have not yet been added as bands) must be added
 * as bands now that they are populated with data and are thus allowed
 * to be added. 
 *
 * If the dataset output bands are fully in memory, they are moved to 
 * a vector from their metadata objects to be passed as a parameter
 * to the GDALRasterWrapper constructor (or not if the bands aren't in
 * memory). This GDALRasterWrapper is then returned.
 *
 * @param GDALRasterWrapper *p_raster
 * @param std::map<int, std::vector<double>>breaks
 * @param bool map
 * @param std::string
 * @param bool largeRaster
 * @param int threads
 * @param std::string tempFolder
 * @param std::map<std::string, std::string> driverOptions
 *
 * p_raster: 
 * 	a pointer to the input raster.
 * breaks: 
 * 	a map of raster band indexes to breaks values.
 * map: 
 * 	a specification of whether to map all of the output
 * 	bands to an additional output band.
 * filename:
 * 	the output filename, or "" if not to write to an 
 * 	output file.
 * largeRaster:
 * 	whether or not the entire raster band should be
 * 	allocated to memory at once.
 * threads:
 * 	the number of threads to process with, only
 * 	used if largeRaster is true.
 * tempFolder:
 * 	the temporary folder to put VRT bands into.
 * driverOptions:
 * 	extra user-defined driver options such as
 * 	compression.
 */
GDALRasterWrapper *breaks(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	std::string filename,
	bool largeRaster,
	int threads,
	std::string tempFolder,
	std::map<std::string, std::string> driverOptions)
{
	GDALAllRegister();

	int height = p_raster->getHeight();
	int width = p_raster->getWidth();
	double *geotransform = p_raster->getGeotransform();
	std::string projection = std::string(p_raster->getDataset()->GetProjectionRef());

	GDALDataset *p_dataset = nullptr;

	size_t bandCount = breaks.size();
	std::vector<std::vector<double>> bandBreaks;
	std::vector<std::string> bandNames = p_raster->getBands();

	std::vector<RasterBandMetaData> dataBands(bandCount);
	std::vector<RasterBandMetaData> stratBands(bandCount + map);
	std::vector<VRTBandDatasetInfo> VRTBandInfo;

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";

	std::mutex dataBandMutex;
	std::mutex stratBandMutex;
	//VRT bands are each their own dataset so they can each have their own mutex
	std::vector<std::mutex> stratBandMutexes(isVRTDataset * (bandCount + map));

	std::string driver;
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

	//allocate, read, and initialize raster data and breaks information
	GDALDataType stratPixelType = GDT_Int8;
	size_t stratPixelSize = 1;
	size_t band = 0;
	for (auto const& [key, val] : breaks) {
		RasterBandMetaData *p_dataBand = &dataBands[band];
		RasterBandMetaData *p_stratBand = &stratBands[band];

		//get and store metadata from input raster band
		GDALRasterBand *p_band = p_raster->getRasterBand(key);
		p_dataBand->p_band = p_band;
		p_dataBand->type = p_raster->getRasterBandType(key);
		p_dataBand->size = p_raster->getRasterBandTypeSize(key);
		p_dataBand->p_buffer = largeRaster ? nullptr : p_raster->getRasterBandBuffer(key);
		p_dataBand->nan = p_band->GetNoDataValue();
		p_dataBand->p_mutex = &dataBandMutex;
		p_band->GetBlockSize(&p_dataBand->xBlockSize, &p_dataBand->yBlockSize);

		//sort and add band breaks vector
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);
		
		//update metadata of new strat raster
		size_t maxStrata = val.size() + 1;
		setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_" + bandNames[key];
		p_stratBand->xBlockSize = map ? dataBands[0].xBlockSize : p_dataBand->xBlockSize;
		p_stratBand->yBlockSize = map ? dataBands[0].yBlockSize : p_dataBand->yBlockSize;
		p_stratBand->p_mutex = &stratBandMutex; //overwritten if VRT dataset

		//update dataset with new band information
		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, *p_stratBand);
		}
		else if (isVRTDataset) {
			createVRTBandDataset(p_dataset, *p_stratBand, tempFolder, std::to_string(key), VRTBandInfo, driverOptions); 
			p_stratBand->p_mutex = &stratBandMutexes[band];
		}
		else { //non-virtual dataset driver
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}

		band++;
	}

	//set multipliers if mapped stratification	
	std::vector<size_t> multipliers(breaks.size(), 1);
	if (map) {
		//determine the stratification band index multipliers of the mapped band
		for (size_t i = 0; i < bandCount - 1; i++) {
			multipliers[i + 1] = multipliers[i] * (bandBreaks[i].size() + 1);
		}

		//update info of new strat raster map band
		RasterBandMetaData *p_stratBand = &stratBands.back();
		size_t maxStrata = multipliers.back() * (bandBreaks.back().size() + 1);
		setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_map";
		p_stratBand->xBlockSize = dataBands[0].xBlockSize;
		p_stratBand->yBlockSize = dataBands[0].yBlockSize;
		p_stratBand->p_mutex = &stratBandMutex; //overwritten if VRT dataset

		//update dataset with new band information
		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, *p_stratBand);			
		}
		else if (isVRTDataset) {
			createVRTBandDataset(
				p_dataset, 
				*p_stratBand, 
				tempFolder, 
				"map", 
				VRTBandInfo, 
				driverOptions
			); 
			p_stratBand->p_mutex = &stratBandMutexes.back();
		}
		else { //non-virtual dataset driver
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}
	}	

	//create full non-virtual dataset now that we have all required band information
	if (!isMEMDataset && !isVRTDataset) {
		bool useTiles = stratBands[0].xBlockSize != width &&
				stratBands[0].yBlockSize != height;

		for (size_t band = 0; band < stratBands.size(); band++) {
			stratBands[band].size = stratPixelSize;
			stratBands[band].type = stratPixelType;
			stratBands[band].p_buffer = !largeRaster ? VSIMalloc3(height, width, stratPixelSize) : nullptr;
		}

		p_dataset = createDataset(
			filename,
			driver,
			width,
			height,
			geotransform,
			projection,
			stratBands.data(),
			stratBands.size(),
			useTiles,
			driverOptions
		);
	}

	//iterate through all pixels and update the stratified raster bands
	if (largeRaster) {
		pybind11::gil_scoped_acquire acquire;
		boost::asio::thread_pool pool(threads);

		if (map) {
			//use the first raster band to determine block size
			int xBlockSize = dataBands[0].xBlockSize;
			int yBlockSize = dataBands[0].yBlockSize;

			int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
			int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;
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
					&dataBands, 
					&stratBands, 
					&bandBreaks,
					&multipliers
				] {
					std::vector<void *> dataBuffers(dataBands.size());
					std::vector<void *> stratBuffers(stratBands.size());
					for (size_t band = 0; band < dataBuffers.size(); band++) {
						dataBuffers[band] = VSIMalloc3(xBlockSize, yBlockSize, dataBands[band].size);
						stratBuffers[band] = VSIMalloc3(xBlockSize, yBlockSize, stratBands[band].size);
					}
					stratBuffers.back() = VSIMalloc3(xBlockSize, yBlockSize, stratBands.back().size);

					for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
						for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
							int xValid, yValid;
							dataBands[0].p_mutex->lock();
							dataBands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
							dataBands[0].p_mutex->unlock();

							//read raster band data into band buffers
							for (size_t band = 0; band < bandCount; band++) {
								rasterBandIO(
									dataBands[band], 
									dataBuffers[band], 
									xBlockSize, 
									yBlockSize, 
									xBlock, 
									yBlock, 
									xValid, 
									yValid, 
									true //read = true
								);		
							}

							//process blocked band data
							for (int y = 0; y < yValid; y++) {
								size_t index = static_cast<size_t>(y * blockSize);
								for (int x = 0; x < xValid; x++) {
									bool mapNan = false;
									size_t mapStrat = 0;
									
									for (size_t band = 0; band < bandCount; band++) {
										processMapPixel(
											index, 
											dataBands[band], 
											dataBuffers[band], 
											stratBands[band], 
											stratBuffers[band], 
											bandBreaks[band],
											multipliers[band],
											mapNan,
											mapStrat
										);
									}
								
									setStrataPixelDependingOnType(
										stratBands.back().type,
										stratBuffers.back(),
										index,
										mapNan,
										mapStrat
									);

									index++;
								}
							}
					
							//write strat band data
							for (size_t band = 0; band <= bandCount; band++) {
								rasterBandIO(
									stratBands[band],
									stratBuffers[band],
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
					}

					for (size_t band = 0; band < dataBuffers.size(); band++) {
						VSIFree(dataBuffers[band]);
						VSIFree(stratBuffers[band]);
					}
					VSIFree(stratBuffers.back());
				});
			}	
		}
		else {	
			for (size_t band = 0; band < bandCount; band++) {
				RasterBandMetaData* p_dataBand = &dataBands[band];
				RasterBandMetaData* p_stratBand = &stratBands[band];

				int xBlockSize = p_dataBand->xBlockSize;
				int yBlockSize = p_dataBand->yBlockSize;
					
				int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
				int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;			
				int chunkSize = yBlocks / threads;
				
				for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
					int yBlockEnd = std::min(yBlocks, yBlockStart + chunkSize);
					std::vector<double> *p_breaks = &bandBreaks[band];
					
					boost::asio::post(pool, [
						xBlockSize, 
						yBlockSize, 
						yBlockStart, 
						yBlockEnd, 
						xBlocks, 
						p_dataBand, 
						p_stratBand, 
						p_breaks
					] {
						void *p_data = VSIMalloc3(xBlockSize, yBlockSize, p_dataBand->size);
						void *p_strat = VSIMalloc3(xBlockSize, yBlockSize, p_stratBand->size);

						int xValid, yValid;
						for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
							for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
								p_dataBand->p_mutex->lock();
								p_dataBand->p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
								p_dataBand->p_mutex->unlock();
								
								//read block into memory
								rasterBandIO(
									*p_dataBand,
									p_data,
									xBlockSize,
									yBlockSize,
									xBlock,
									yBlock,
									xValid,
									yValid,
									true //read = true
								);

								//process block
								for (int y = 0; y < yValid; y++) {
									index = static_cast<size_t>(y * xblockSize);
									for (int x = 0; x < xValid; x++) {
										processPixel(index, p_data, p_dataBand, p_strat, p_stratBand, *p_breaks);
										index++;
									}
								}
								
								//write resulting stratifications to disk
								rasterBandIO(
									*p_stratBand,
									p_strat,
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
						VSIFree(p_data);
						VSIFree(p_strat);
					});
				}
			}
		}
		pool.join();
		pybind11::gil_scoped_release release;
	}
	else {
		size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) * static_cast<size_t>(p_raster->getHeight());
		if (map) {
			for (size_t index = 0; index < pixelCount; index++) {
				bool mapNan = false;
				size_t mapStrat = 0;
									
				for (size_t band = 0; band < bandCount; band++) {
					processMapPixel(
						index, 
						dataBands[band], 
						dataBands[band].p_buffer, 
						stratBands[band], 
						stratBands[band].p_buffer, 
						bandBreaks[band],
						multipliers[band],
						mapNan,
						mapStrat
					);
				}

				setStrataPixelDependingOnType(
					stratBands.back().type,
					stratBands.back().p_buffer,
					index,
					mapNan,
					mapStrat
				);
			}
		}
		else {
			for (size_t band = 0; band < bandCount; band++) {
				for (size_t i = 0; i < pixelCount; i++) {
					processPixel(
						i,
					       	dataBands[band].p_buffer,	
						&dataBands[band], 
						stratBands[band].p_buffer,
						&stratBands[band],
						bandBreaks[band]
					);
				}
			}
		}

		//if non-virtual dataset then write in-memory output bands to disk
		if (!isVRTDataset && !isMEMDataset) {
			CPLErr err;
			for (size_t band = 0; band < stratBands.size(); band++) {
				err = stratBands[band].p_band->RasterIO(
					GF_Write,
					0,
					0,
					width,
					height,
					stratBands[band].p_buffer,
					width,
					height, 
					stratBands[band].type,
					0,
					0
				);
				if (err) {
					throw std::runtime_error("error writing band to file.");
				}
			}
		}
	}

	//close and add all of the VRT sub datasets as bands
	if (isVRTDataset) {
		for (size_t band = 0; band < VRTBandInfo.size(); band++) {
			GDALClose(VRTBandInfo[band].p_dataset);
			addBandToVRTDataset(p_dataset, stratBands[band], VRTBandInfo[band]);
		}
	}

	//if the bands are in memory, populate a vector of just bands to use in GDALRasterWrapper creation
	std::vector<void *> buffers(stratBands.size());
	if (!largeRaster) {
		for (size_t band = 0; band < stratBands.size(); band++) {
			buffers[band] = stratBands[band].p_buffer;
		}
	}

	return largeRaster ?
	       	new GDALRasterWrapper(p_dataset) :
		new GDALRasterWrapper(p_dataset, buffers);

}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaks);
}
