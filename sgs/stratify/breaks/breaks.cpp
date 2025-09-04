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

#define GIGABYTE 1073741824

/**
 *
 */
inline void processMapPixel(
	size_t index,
	std::vector<RasterBandMetaData>& dataBands,
	std::vector<void *>& dataBuffers,
	std::vector<RasterBandMetaData>& stratBands,
	std::vector<void *>& stratBuffers,
	std::vector<size_t>& multipliers,
	std::vector<std::vector<double>>& bandBreaks)
{
	size_t mappedStrat = 0;
	bool mapNan = false;
	for (size_t band = 0; band < dataBands.size(); band++) {
		double val = getPixelValueDependingOnType<double>(
			dataBands[band].type,
			dataBuffers[band],
			index
		);
		bool isNan = std::isnan(val) || (double)val == dataBands[band].nan;
		mapNan |= isNan;
		
		//calculate strata value if not nan
		size_t strat = 0;
		if (!isNan) {
			std::vector<double> curBandBreaks = bandBreaks[band];
			auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
			strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
		}

		setStrataPixelDependingOnType(
			stratBands[band].type,
			stratBuffers[band],
			index,
			isNan,
			strat
		);

		//adjust mappedStrat as required
		if (!mapNan) {
			mappedStrat += strat * multipliers[band];
		}
	}
		
	//assign mapped value
	setStrataPixelDependingOnType(
		stratBands.back().type,
		stratBuffers.back(),
		index,
		mapNan,
		mappedStrat
	);
}	

inline void
processMapPixel(
	size_t index,
	std::vector<RasterBandMetaData>& dataBands,
	std::vector<RasterBandMetaData>& stratBands,
	std::vector<size_t>& multipliers,
	std::vector<std::vector<double>>& bandBreaks)
{
	std::vector<void *> dataBuffers;
	std::vector<void *> stratBuffers;
	for (const RasterBandMetaData& dataBand : dataBands) {
		dataBuffers.push_back(dataBand.p_buffer);
	}
	for (const RasterBandMetaData& stratBand : stratBands) {
		stratBuffers.push_back(stratBand.p_buffer);
	}
	
	processMapPixel(
		index, 
		dataBands, 
		dataBuffers, 
		stratBands, 
		stratBuffers, 
		multipliers, 
		bandBreaks
	);
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
	bool isNan = std::isnan(val) || (double)val == p_dataBand->nan;

	//calculate strata value if not nan
	size_t strat = 0;
	if (!isNan) {
		std::vector<double> curBandBreaks = bandBreaks;
		auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
		strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
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
inline void
processPixel(
	size_t index,
	RasterBandMetaData *p_dataBand,
	RasterBandMetaData *p_stratBand,
	std::vector<double>& bandBreaks
) {
	processPixel(
		index, 
		p_dataBand->p_buffer, 
		p_dataBand, 
		p_stratBand->p_buffer, 
		p_stratBand, 
		bandBreaks
	);
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
			addBandToMEMDataset(p_dataset, p_stratBand);
		}
		else if (isVRTDataset) {
			std::filesystem::path tmpPath = tempFolder;
			std::filesystem::path tmpName = "strat_breaks_" + std::to_string(key) + ".tif";
			tmpPath = tmpPath / tmpName;

			VRTBandDatasetInfo info;		
			info.filename = tmpPath.string();
			createVRTSubDataset(p_dataset, p_stratBand, info);
			VRTBandInfo.push_back(info);
		}
		else { //dataset driver is the one which corresponds to the filename extension
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
			addBandToMEMDataset(p_dataset, p_stratBand);
		}
		else if (isVRTDataset) {
			std::filesystem::path tmpPath = tempFolder;
			std::filesystem::path tmpName = "strat_breaks_map.tif";
			tmpPath = tmpPath / tmpName;

			VRTBandDatasetInfo info;
			info.filename = tmpPath.string();
			createVRTSubDataset(p_dataset, p_stratBand, info);
			VRTBandInfo.push_back(info);
		}
		else { //dataset driver is the one which corresponds to the filename extension
			if (stratPixelSize < p_stratBand->size) {
				stratPixelType = p_stratBand->type;
				stratPixelSize = p_stratBand->size;
			}
		}
	}	

	//now we can create the dataset if the type is not MEM or VRT because we know the types of each band
	if (!isMEMDataset && !isVRTDataset) {
		//tiles must not be scanlines, as trying to set block size when they represent scanlines 
		//may result in GDAL errors due to the tile array being too large. 
		bool useTiles = largeRaster && 
				stratBands[0].xBlockSize != width && 
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
			useTiles
		);
	}

	//step 3: iterate through indices and update the stratified raster bands
	if (largeRaster) {
		pybind11::gil_scoped_acquire acquire;
		unsigned int threads = std::min(static_cast<unsigned int>(8), std::thread::hardware_concurrency());
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
							CPLErr err;
							for (size_t band = 0; band < bandCount; band++) {
								bool useReadBlock = xBlockSize == dataBands[0].xBlockSize &&
										    yBlockSize == dataBands[0].yBlockSize;
								
								dataBands[band].mutex.lock();
								if (useReadBlock) {
									err = dataBands[band].p_band->ReadBlock(
										xBlock, 
										yBlock, 
										dataBuffers[band]
									);
								}
								else {
									err = dataBands[band].p_band->RasterIO(
										GF_Read,
										xBlock * xBlockSize,
										yBlock * yBlockSize,
										xValid,
										yValid,
										dataBuffers[band],
										xBlockSize,
										yBlockSize,
										dataBands[band].type,
										0,
										0
									);
								}
								dataBands[band].mutex.unlock();
								if (err) {
									throw std::runtime_error("unable to read block from raster");
								}
							}

							//process blocked band data
							for (int y = 0; y < yValid; y++) {
								for (int x = 0; x < xValid; x++) {
									size_t index = static_cast<size_t>(x + y * xBlockSize);
									processMapPixel(index, dataBands, dataBuffers, stratBands, stratBuffers, multipliers, bandBreaks);
								}
							}
					
							//write strat band data
							for (size_t band = 0; band <= bandCount; band++) {
								bool useWriteBlock = xBlockSize == stratBands[0].xBlockSize &&
										     yBlockSize == stratBands[0].yBlockSize;
								
								stratBands[band].mutex.lock();
								if (useWriteBlock) {
									err = stratBands[band].p_band->WriteBlock(
										xBlock, 
										yBlock, 
										stratBuffers[band]
									);
								}
								else {
									err = stratBands[band].p_band->RasterIO(
										GF_Write,
										xBlock * xBlockSize,
										yBlock * yBlockSize,
										xValid,
										yValid,
										stratBuffers[band],
										xBlockSize,
										yBlockSize,
										stratBands[band].type,
										0,
										0
									);
								}
								stratBands[band].mutex.unlock();
								if (err) {
									throw std::runtime_error("error writing block of data to band.");
								}
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
					
				if (xBlockSize != p_stratBand->xBlockSize || 
				    yBlockSize != p_stratBand->yBlockSize) {
					throw std::runtime_error("block size error -- should not be here!!!");
				}

				int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
				int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;			
				int chunkSize = yBlocks / threads;
				
				for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
					int yBlockEnd = std::min(yBlocks, yBlockStart + chunkSize);
					std::vector<double> breakVals = bandBreaks[band];
					
					boost::asio::post(pool, [
						xBlockSize, 
						yBlockSize, 
						yBlockStart, 
						yBlockEnd, 
						xBlocks, 
						p_dataBand, 
						p_stratBand, 
						&breakVals
					] {
						bool useWriteBlock = p_stratBand->xBlockSize == xBlockSize &&
										p_stratBand->yBlockSize == yBlockSize;
						
						CPLErr err;
						void *p_data = VSIMalloc3(xBlockSize, yBlockSize, p_dataBand->size);
						void *p_strat = VSIMalloc3(xBlockSize, yBlockSize, p_stratBand->size);
						int xValid, yValid;
						for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
							for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
								p_dataBand->mutex.lock();
								err = p_dataBand->p_band->ReadBlock(xBlock, yBlock, p_data);
								p_dataBand->p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
								p_dataBand->mutex.unlock();
								if (err) {
									throw std::runtime_error("error reading block from band.");
								}

								for (int y = 0; y < yValid; y++) {
									for (int x = 0; x < xValid; x++) {
										size_t index = static_cast<size_t>(x + y * xBlockSize);
										processPixel(index, p_data, p_dataBand, p_strat, p_stratBand, breakVals);
									}
								}

								p_stratBand->mutex.lock();
								if (useWriteBlock) {
									err = p_stratBand->p_band->WriteBlock(xBlock, yBlock, p_strat);
								}
								else {
									err = p_stratBand->p_band->RasterIO(
										GF_Write,
										xBlock * xBlockSize,
										yBlock * yBlockSize,
										xValid,
										yValid,
										p_strat,
										xBlockSize,
										yBlockSize,
										p_stratBand->type,
										0,
										0
									);
								}
								p_stratBand->mutex.unlock();
								if (err) {
									throw std::runtime_error("error writing block to band.");
								}
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
			for (size_t i = 0; i < pixelCount; i++) {
				processMapPixel(i, dataBands, stratBands, multipliers, bandBreaks);		
			}
		}
		else {
			for (size_t band = 0; band < bandCount; band++) {
				for (size_t i = 0; i < pixelCount; i++) {
					processPixel(
						i, 
						&dataBands[band], 
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
		
	return p_stratRaster;
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaks);
}
