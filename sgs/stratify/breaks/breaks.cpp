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
 *
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
	double val = getPixelValueDependingOnType<double>(
		dataBand.type,
		p_dataBuffer,
		index
	);
	bool isNan = std::isnan(val) || (double)val == dataBand.nan;
	mapNan |= isNan;
		
	//calculate strata value if not nan
	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(bandBreaks.begin(), bandBreaks.end(), val);
		strat = (it == bandBreaks.end()) ? bandBreaks.size() : std::distance(bandBreaks.begin(), it);
	}

	setStrataPixelDependingOnType(
		stratBand.type,
		p_stratBuffer,
		index,
		isNan,
		strat
	);

	//adjust mappedStrat as required
	if (!mapNan) {
		mapStrat += strat * multiplier;
	}
}	

/**
 *
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
	double val = getPixelValueDependingOnType<double>(
		p_dataBand->type,
		p_data,
		index
	);

	bool isNan = std::isnan(val) || val == p_dataBand->nan;

	//calculate strata value if not nan
	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(bandBreaks.begin(), bandBreaks.end(), val);
		strat = (it == bandBreaks.end()) ? bandBreaks.size() : std::distance(bandBreaks.begin(), it);
	}

	setStrataPixelDependingOnType(
		p_stratBand->type,
		p_strat,
		index,
		isNan,
		strat
	);
}

/**
 *
 */
template <typename T, typename U>
inline void 
processBand(
	GDALRasterWrapper *p_raster,
	RasterBandMetaData& dataBand,
	RasterBandMetaData& stratBand,
	boost::asio::thread_pool& pool,
	std::vector<double>& bandBreaks,
	size_t threads)
{
	int xBlockSize = dataBand.xBlockSize;
	int yBlockSize = dataBand.yBlockSize;

	int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
	int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;
	int chunkSize = yBlocks / threads;
	
	for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
		int yBlockEnd = std::min(yBlocks, yBlockStart + chunkSize);
		
		boost::asio::post(pool, [
			xBlockSize, 
			yBlockSize, 
			yBlockStart, 
			yBlockEnd, 
			xBlocks, 
			&dataBand, 
			&stratBand,
			&bandBreaks 
		] {
			//these will be freed by the thread which uses them
			T *p_data = reinterpret_cast<T *>(VSIMalloc3(xBlockSize, yBlockSize, sizeof(T)));
			U *p_strat = reinterpret_cast<U *>(VSIMalloc3(xBlockSize, yBlockSize, sizeof(U)));

			T dataNan = static_cast<T>(dataBand.nan);
			U stratNan = static_cast<U>(stratBand.nan);
			std::vector<T> breaks = std::vector<T>(bandBreaks.size());
			for (size_t i = 0; i < bandBreaks.size(); i++) {
				breaks[i] = static_cast<T>(bandBreaks[i]);
			}

			std::chrono::microseconds read = std::chrono::microseconds(0);
			std::chrono::microseconds process = std::chrono::microseconds(0);
			std::chrono::microseconds write = std::chrono::microseconds(0);

			int xValid, yValid;
			for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
				for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
					dataBand.mutex.lock();
					dataBand.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
					dataBand.mutex.unlock();

					const auto read1 = std::chrono::high_resolution_clock::now();
					rasterBandIO(
						dataBand,
						p_data,
						xBlockSize,
						yBlockSize,
						xBlock,
						yBlock,
						xValid,
						yValid,
						true //read = true
					);
					const auto read2 = std::chrono::high_resolution_clock::now();
					read += std::chrono::duration_cast<std::chrono::microseconds>(read2 - read1);

					const auto process1 = std::chrono::high_resolution_clock::now();
					for (int y = 0; y < yValid; y++) {
						for (int x = 0; x < xValid; x++) {
							int index = x + y * xBlockSize;
							T val = p_data[index];
							bool isNan = val == dataNan;
							if (!isNan) {
								auto it = std::lower_bound(breaks.begin(), breaks.end(), val);
								size_t strat = (it == breaks.end()) ? breaks.size() : std::distance(breaks.begin(), it);
								p_strat[index] = strat;
							}
							else {
								p_strat[index] = stratNan;
							}
						}
					}
					const auto process2 = std::chrono::high_resolution_clock::now();
					process += std::chrono::duration_cast<std::chrono::microseconds>(process2 - process1);

					const auto write1 = std::chrono::high_resolution_clock::now();					
					rasterBandIO(
						stratBand,
						p_strat,
						xBlockSize,
						yBlockSize,
						xBlock,
						yBlock,
						xValid,
						yValid,
						false //read = false
					);
					const auto write2 = std::chrono::high_resolution_clock::now();
					write += std::chrono::duration_cast<std::chrono::microseconds>(write2 - write1);
				}
			}
			VSIFree(p_data);
			VSIFree(p_strat);

			std::cout << "total read time: " << std::chrono::duration_cast<std::chrono::seconds>(read).count() << std::endl;
			std::cout << "total process time: " << std::chrono::duration_cast<std::chrono::seconds>(process).count() << std::endl;
			std::cout << "total write time: " << std::chrono::duration_cast<std::chrono::seconds>(write).count() << std::endl;
		});
	}
}

/**
 * This function stratifies a given raster using user-defined breaks.
 * The breaks are provided as a vector of doubles mapped to a band index.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map variable to combine the stratification of
 * the multiple raster bands.
 *
 * The required raster bands are first acquired from the datset, and
 * checked to ensure that the number of breaks fits without overflowing.
 *
 * An additional band is added to the output raster if map is specified,
 * and multipliers are determined, so that every unique combination
 * of stratifications of the normal raster bands corresponds to a 
 * single stratification of the mapped raster band. For example,
 * if there were 3 normal raster bands each with 5 possible stratifications,
 * there would be 5^3 or 125 possible mapped stratificaitons.
 *
 * The raster bands are then iterated thorough and stratifications are determined.
 * a new dataset object is created using the stratifcation raster, and it is
 * written to disk if a filename is given.
 *
 * NOTE: stratifications have a data type of 'float' to accomodate nan values
 *
 * @param GDALRasterWrapper * a pointer to the raster image to stratify
 * @param std::map<int, std::vector<double>> band mapped to user-defiend breaks
 * @param bool map whether to add a mapped stratification
 * @param std::string filename the filename to write to (if desired)
 */
GDALRasterWrapper *breaks(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	std::string filename,
	bool largeRaster,
	std::string tempFolder)
{
	const auto setupStart = std::chrono::high_resolution_clock::now();
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

	//step 1: allocate, read, and initialize raster data and breaks information
	GDALDataType stratPixelType = GDT_Int8;
	size_t stratPixelSize = 1;
	size_t band = 0;
	for (auto const& [key, val] : breaks) {
		RasterBandMetaData *p_dataBand = &dataBands[band];
		RasterBandMetaData *p_stratBand = &stratBands[band];

		GDALRasterBand *p_band = p_raster->getRasterBand(key);
		p_dataBand->p_band = p_band;
		p_dataBand->type = p_raster->getRasterBandType(key);
		p_dataBand->size = p_raster->getRasterBandTypeSize(key);
		p_dataBand->p_buffer = largeRaster ? nullptr : p_raster->getRasterBandBuffer(key);
		p_dataBand->nan = p_band->GetNoDataValue();
		p_band->GetBlockSize(&p_dataBand->xBlockSize, &p_dataBand->yBlockSize);

		//sort and add band breaks vector
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);
		
		//update info of new strat raster
		size_t maxStrata = val.size() + 1;
		setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_" + bandNames[key];
		p_stratBand->xBlockSize = map ? dataBands[0].xBlockSize : p_dataBand->xBlockSize;
		p_stratBand->yBlockSize = map ? dataBands[0].yBlockSize : p_dataBand->yBlockSize;
		
		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, *p_stratBand);			
		}
		else if (isVRTDataset) {
			addBandToVRTDataset(p_dataset, *p_stratBand, tempFolder, std::to_string(key), VRTBandInfo); 
		}
		else { //non-virtual dataset driver
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}

		band++;
	}

	//step 2: set multipliers and check max size if mapped stratification	
	std::vector<size_t> multipliers(breaks.size(), 1);
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
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

		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, *p_stratBand);			
		}
		else if (isVRTDataset) {
			addBandToVRTDataset(p_dataset, *p_stratBand, tempFolder, "map", VRTBandInfo); 
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
		p_dataset = adjustBandsAndCreateDataset(
			filename,
			driver,
			width,
			height,
			geotransform,
			projection,
			stratBands,
			stratPixelSize,
			stratPixelType,
			largeRaster	
		);		
	}
	std::cout << "setup takes: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - setupStart).count() << std::endl;

	const auto iterateStart = std::chrono::high_resolution_clock::now();
	//step 3: iterate through indices and update the stratified raster bands
	if (largeRaster) {
		pybind11::gil_scoped_acquire acquire;
		unsigned int threads = std::min(static_cast<unsigned int>(1), std::thread::hardware_concurrency());
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
							dataBands[0].mutex.lock();
							dataBands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
							dataBands[0].mutex.unlock();

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
								for (int x = 0; x < xValid; x++) {
									size_t index = static_cast<size_t>(x + y * xBlockSize);
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
			if (threads > 1) {
				throw std::runtime_error("chrono runtime tests won't work because of interleaving threads!!");
			}

			for (size_t band = 0; band < bandCount; band++) {
				int dataBandType = static_cast<int>(dataBands[band].type);
				int stratBandType = static_cast<int>(stratBands[band].type);
				switch((dataBandType << 4) + stratBandType) {
					case (16 * 14 + 14): {processBand<int8_t, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 14 + 3): {processBand<int8_t, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 14 + 5): {processBand<int8_t, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 2 + 14): {processBand<uint16_t, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 2 + 3): {processBand<uint16_t, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 2 + 5): {processBand<uint16_t, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 3 + 14): {processBand<int16_t, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 3 + 3): {processBand<int16_t, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 3 + 5): {processBand<int16_t, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 4 + 14): {processBand<uint32_t, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 4 + 3): {processBand<uint32_t, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 4 + 5): {processBand<uint32_t, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 5 + 14): {processBand<int32_t, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 5 + 3): {processBand<int32_t, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 5 + 5): {processBand<int32_t, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 6 + 14): {processBand<float, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 6 + 3): {processBand<float, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 6 + 5): {processBand<float, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 7 + 14): {processBand<double, int8_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 7 + 3): {processBand<double, int16_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					case (16 * 7 + 5): {processBand<double, int32_t>(p_raster, dataBands[band], stratBands[band], pool, bandBreaks[band], threads); break;}
					default: {throw std::runtime_error("raster band data type not supported!");}
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
	std::cout << "iterate takes: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - iterateStart).count() << std::endl;

	const auto closeStart = std::chrono::high_resolution_clock::now();	
	//free allocated band data
	if (map && largeRaster) {
		for (size_t band = 0; band < stratBands.size(); band++) {
			VSIFree(stratBands[band].p_buffer);
		}
	}

	//step 4: close all of the VRT sub datasets
	if (isVRTDataset) {
		for (size_t band = 0; band < VRTBandInfo.size(); band++) {
			GDALClose(VRTBandInfo[band].p_dataset);
			addBandToVRTDataset(p_dataset, stratBands[band], VRTBandInfo[band]);
		}
	}

	//step 5: if the bands are in memory, populate a vector of just bands to use in GDALRasterWrapper creation
	std::vector<void *> buffers(stratBands.size());
	if (!largeRaster) {
		for (size_t band = 0; band < stratBands.size(); band++) {
			buffers[band] = stratBands[band].p_buffer;
		}
	}

	//step 5: create GDALRasterWrapper object from bands
	//this dynamically-allocated object will be cleaned up by python (TODO I hope...)
	GDALRasterWrapper *p_stratRaster = largeRaster ?
	       	new GDALRasterWrapper(p_dataset) :
		new GDALRasterWrapper(p_dataset, buffers);
		

	std::cout << "close takes: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - closeStart).count() << std::endl;
	return p_stratRaster;
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaks);
}
