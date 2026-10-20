/******************************************************************************
 *
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using breaks
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

/**
 * @defgroup breaks breaks
 * @ingroup stratify
 */

#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "utils/raster.h"
#include "utils/helper.h"

namespace sgs {
namespace breaks {

{	
	for (int y = 0; y < yValid; y++) {
		size_t index = static_cast<size_t>(y * xBlockSize);
		for (int x = 0; x < xValid; x++) {
			bool mapNan = false;
			T mapStrat = 0;
			
			for (size_t band = 0; band < bandCount; band++) {
				D val = dataBuffers[band][index];
				bool isNan = std::isnan(val) || val == bandsnan[band];
				mapNan |= isNan;
	
				if (!isNan) {
					auto it = std::lower_bound(bandBreaks[band].begin(), bandBreaks[band].end(), val);
					size_t strat = (it == bandBreaks[band].end()) ? bandBreaks[band].size() : std::distance(bandBreaks[band].begin(), it);
					stratBuffers[band][index] = static_cast<S>(strat);
				}
				else {
					stratBuffers[band][index] = stratNan;
				}
		
				if (!mapNan) {
					mapStrat += strat * multipliers[band];
				}
			}
			
			stratBuffers.back()[index] = mapNan ? stratNan : mapStrat;
			index++;
		}
	}
}

template <typename D, typename S>
void
stratifyMapped(
	std::vector<RasterBandMetaData>& dataBands,
	std::vector<RasterBandMetaData>& stratBands,
	GDALDataType dataBandType,
	GDALDataType stratBandType,
	int height,
	int width,
	std::vector<std::vector<D>>& bandBreaks,
	std::vector<size_t>& multipliers,
	int threads) 
{
	pybind11::gil_scoped_acquire acquire;
	boost::asio::thread_pool pool(threads);

	std::vector<S> sMultipliers(multipliers.size());
	for (size_t i = 0; i < multipliers.size(); i++) {
		sMultipliers[i] = static_cast<S>(multipliers[i]);
	}

	std::vector<D> bandsNan(dataBands.size());
	for (size_t i = 0; i < dataBands.size(); i++) {
		bandsNan[i] = static_cast<D>(dataBands[i].nan);
	}

	S stratNan = -1;

	//use the first raster band to determine block size
	int xBlockSize = dataBands[0].xBlockSize;
	int yBlockSize = dataBands[0].yBlockSize;

	int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
	int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;
	int chunkSize = yBlocks / threads;

	for (int yBlockStart = 0; yBlockStart < yBlocks; yBlockStart += chunkSize) {
		int yBlockEnd = std::min(yBlockStart + chunkSize, yBlocks);
	
		boost::asio::post(pool, [
			height,
			width,
			xBlockSize, 
			xBlocks, 
			yBlockSize, 
			yBlockStart, 
			yBlockEnd, 
			&dataBands, 
			&stratBands, 
			&bandBreaks,
			&sMultipliers
			&bandsNan,
			stratNan
		] {
			
			if (allDataInMemory) {
				for (int y = yBlockStart * yBlockSize; y < yBlockEnd * yBlockSize; y++) {
					size_t i = static_cast<size_t>(y * xBlockSize);
					for (int x = 0; x < width; x++) {
						bool mapNan = false;
						T mapStrat = 0;

						for (size_t b = 0; b < dataBands.size(); b++) {
							D val = reinterpret_cast<D *>(dataBuffers[b])[i];
							bool isNan = std::isnan(val) || val == bandsNan[b];
							mapNan |= isNan;

							if (!isNan) {
								auto it = std::lower_bound(
									bandBreaks[b].begin(), 
									bandBreaks[b].end(), 
									val);
								size_t strat = (it == bandBreaks[b].end()) ? 
										bandBreaks[b].size() : 
										std::distance(bandBreaks[b].begin(), it);
								stratBuffers[b][i] = static_cast<S>(strat);
							}
							else {
								stratBuffers[b][i] = stratNan;
							}

							if (!mapNan) {
								mapStrat += strat * sMultipliers[band];
							}
						}

						stratBuffers.back()[index] = mapNan ? stratNan : mapStrat;
						i++;
					}
				}
			}
			else { //!allDataInMemory
				std::vector<D *> dataBuffers(dataBands.size());
				std::vector<S *> stratBuffers(stratBands.size());

				for (size_t band = 0; band < dataBuffers.size(); band++) {
					dataBuffers[band] = static_cast<D *>(
						VSIMalloc3(xBlockSize, yBlockSize, dataPixelSize));
					stratBuffers[band] = static_cast<S *>(
						VSIMalloc3(xBlockSize, yBlockSize, stratBands[band].size));
				}	
				stratBuffers.back() = static_cast<S *>(
					VSIMalloc3(xBlockSize, yBlockSize, stratBands.back().size));
				
				//process the raster in blocks, 
				for (int yBlock = yBlockStart; yBlock < yBlockEnd; yBlock++) {
					for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
						int xValid, yValid;
						dataBands[0].p_mutex->lock();
						dataBands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
						dataBands[0].p_mutex->unlock();

						//read raster band data into band buffers
						for (size_t band = 0; band < bandCount; band++) {
							helper::rasterBandIO(
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

						for (int y = 0; y < yValid; y++) {
							size_t i = static_cast<size_t>(y * xBlockSize);
							for (int x = 0; x < xValid; x++) {
								bool mapNan = false;
								T mapStrat = 0;

								for (size_t b = 0; b < dataBands.size(); b++) {
									D val = reinterpret_cast<D *>(dataBuffers[b])[i];
									bool isNan = std::isnan(val) || val == bandsNan[b];
									mapNan |= isNan;

									if (!isNan) {
										auto it = std::lower_bound(
											bandBreaks[b].begin(), 
											bandBreaks[b].end(), 
											val);
										size_t strat = (it == bandBreaks[b].end()) ? 
												bandBreaks[b].size() : 
												std::distance(bandBreaks[b].begin(), it);
										stratBuffers[b][i] = static_cast<S>(strat);
									}
									else {
										stratBuffers[b][i] = stratNan;
									}

									if (!mapNan) {
										mapStrat += strat * sMultipliers[band];
									}
								}

								stratBuffers.back()[index] = mapNan ? stratNan : mapStrat;
								i++;
							}
						}
	
						//write strat band data
						for (size_t band = 0; band <= bandCount; band++) {
							helper::rasterBandIO(
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
			}

			if (!allDataInMemory) {
				for (size_t band = 0; band < dataBuffers.size(); band++) {
					VSIFree(dataBuffers[band]);
					VSIFree(stratBuffers[band]);
				}
				VSIFree(stratBuffers.back());
			}
		});
	}	

	pool.join();
	pybind11::gil_scoped_release release;
}

//TODO UPDATE TO BLOCK
/**
 * @ingroup breaks
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
 *
 * @param size_t index
 * @param void *p_data
 * @param RasterBandMetaData *p_dataBand
 * @param void *p_strat
 * @param RasterBandMetaData *p_stratBand
 * @param std::vector<double>& bandBreaks
 */
inline void
processBlock(
	int xValid,
	int yValid,
	int xBlockSize,
	void *p_dataBuffer,
	helper::RasterBandMetaData& dataBand,
	void *p_stratBuffer,
	helper::RasterBandMetaData& stratBand,
	std::vector<double>& bandBreaks)
{
	for (int y = 0; y < yValid; y++) {
		size_t index = static_cast<size_t>(y * xBlockSize);
		for (int x = 0; x < xValid; x++) {
			double val = helper::getPixelValueDependingOnType<double>(dataBand.type, p_dataBuffer, index);
			bool isNan = std::isnan(val) || val == dataBand.nan;
		
			size_t strat = 0;
			if (!isNan) {
				auto it = std::lower_bound(bandBreaks.begin(), bandBreaks.end(), val);
				strat = (it == bandBreaks.end()) ? bandBreaks.size() : std::distance(bandBreaks.begin(), it);
			}
		
			helper::setStrataPixelDependingOnType(stratBand.type, p_stratBuffer, index, isNan, strat);
			index++;
		}
	}
	
	
}

template <typename D, typename S>
void
stratify() {
	pybind11::gil_scoped_acquire acquire;
	boost::asio::thread_pool pool(threads);

	helper::RasterBandMetaData* p_dataBand = &dataBands[band];
	helper::RasterBandMetaData* p_stratBand = &stratBands[band];

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
					helper::rasterBandIO(
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

					//PROCESS BLOCK
					
					//write resulting stratifications to disk
					helper::rasterBandIO(
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

	pool.join();
	pybind11::gil_scoped_release release;
}

/**
 * @ingroup breaks
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
 * @param std::string filename
 * @param bool largeRaster
 * @param int threads
 * @param std::string tempFolder
 * @param std::map<std::string, std::string> driverOptions
 * @returns GDALRasterWrapper *
 */
raster::GDALRasterWrapper *breaks(
	raster::GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> breaks,
	bool map,
	std::string filename,
	int threads,
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

	std::vector<helper::RasterBandMetaData> dataBands(bandCount);
	std::vector<helper::RasterBandMetaData> stratBands(bandCount + map);
	std::vector<helper::VRTBandDatasetInfo> VRTBandInfo;

	std::mutex dataBandMutex;
	std::mutex stratBandMutex;

	bool onDisk = filename != "";

	GDALDataType stratPixelType = GDT_Int8;
	size_t stratPixelSize = sizeof(uint8_t);
	
	GDALDataType dataPixelType = GDT_Float32;
	size_t dataPixelSize = sizeof(float);
	
	//allocate, read, and initialize raster data and breaks information
	size_t i = 0;
	for (auto const& [key, val] : breaks) {
		//get and store metadata from input raster band
		p_raster->fillRasterBandMetaData(key, dataBands[i]);
		dataBands[i].p_mutex = &dataBandMutex;

		if (dataBands[i].type == GDT_Float64 && dataPixelType == GDT_Float32) {
			dataPixelType = GDT_Float32;
			dataPixelSize = sizeof(double);
		}

		//sort and add band breaks vector
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);
		
		//update metadata of new strat raster
		size_t maxStrata = val.size() + 1;
		helper::setStratBandTypeAndSize(maxStrata, &stratBands[i].type, &stratBands[i].size);
		stratBands[i].name = "strat_" + bandNames[key];
		stratBands[i].xBlockSize = map ? dataBands[0].xBlockSize : dataBands[i].xBlockSize;
		stratBands[i].yBlockSize = map ? dataBands[0].yBlockSize : dataBands[i].yBlockSize;
		stratBands[i].p_mutex = &stratBandMutex; //overwritten if VRT dataset

		//if writing to file on disk, ensure type and size are the largest any single band requires
		if (stratPixelSize < stratBands[i].size) {
			stratPixelSize = stratBands[i].size;
			stratPixelType = stratBands[i].type;
		}

		i++;
	}

	//set multipliers if mapped stratification	
	std::vector<size_t> multipliers(breaks.size(), 1);
	if (map) {
		//determine the stratification band index multipliers of the mapped band
		for (size_t i = 0; i < bandCount - 1; i++) {
			multipliers[i + 1] = multipliers[i] * (bandBreaks[i].size() + 1);
		}

		//update info of new strat raster map band
		size_t index = stratBands.size() - 1;
		stratBands[index].name = "strat_map";
		stratBands[index].xBlockSize = dataBands[0].xBlockSize;
		stratBands[index].yBlockSize = dataBands[0].yBlockSize;
		stratBands[index].p_mutex = &stratBandMutex; //overwritten if VRT dataset
		size_t maxStrata = multipliers.back() * (bandBreaks.back().size() + 1);
		helper::setStratBandTypeAndSize(maxStrata, &stratBands[index].type, &stratBands[index].size);
		
		//type and size must be the same for all bands in an output file because of
		//restrictions on GeoTiff files which aren't restricted in virtual types
		if (stratPixelSize < stratBands[index].size) {
			stratPixelSize = stratBands[index].size;
			stratPixelType = stratBands[index].type;
		}
	}	

	//if on disk (geotif), all bands of the raster MUST have the same type
	//which is not a requirement of virtual file types.
	if (onDisk || map) {
		for (size_t i = 0; i < stratBands.size(); i++) {
			stratBands[i].size = stratPixelSize;
			stratBands[i].type = stratPixelType;
		}
	}

	std::vector<void *> allocated;
	bool allocatedWholeBands = onDisk ? false : 
		helper::wholeBandsAllocated(height, width, dataPixelType, dataPixelSize,
					    dataBands, stratBands, allocated);
	
	bool isMEMDataset = allocatedWholeBands && !onDisk;
	bool isVRTDataset = !allocatedWholeBands && !onDisk;
	
	//construct mutexes for each VRT band dataset if required
	std::vector<std::mutex> stratBandMutexes(isVRTDataset * (bandCount + map));

	if (isMEMDataset) {
		p_dataset = helper::createVirtualDataset("MEM", width, height, geotransform, projection);

		for (helper::RasterBandMetaData& band : stratBands) {
			helper::addBandToMEMDataset(p_dataset, band);
		}
	}
	else if (isVRTDataset) {
		p_dataset = helper::createVirtualDataset("VRT", width, height, geotransform, projection);

		for (size_t i = 0; i < stratBands.size(); i++) {
			helper::createVRTBandDataset(p_dataset, stratBands[i], VRTBandInfo, driverOptions);
			stratBands[i].p_mutex = &stratBandMutexes[i];
		}
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();
		std::string driver;

		if (extension == ".tif") {
			driver = "Gtiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now.");
		}

		int xBlockSize = stratBands[0].xBlockSize;
		int yBlockSize = stratBands[0].yBlockSize;
		bool useTiles = xBlockSize != width &&
				yBlockSize != height;

		p_dataset = helper::createDataset(
			filename,
			driver,
			width,
			height,
			geotransform,
			projection,
			stratBands,
			useTiles,
			driverOptions
		);
	}

	//get pixel types as uint32_t to pass combination to case switch statement
	uint32_t dataBandType = static_cast<uint32_t>(dataPixelType);
	uint32_t stratBandType = static_cast<uint32_t>(stratPixelType);

	//get gdal pixel types as uint32_t to pass combinations to case switch statement
	uint32_t gdalf32 = static_cast<uint32_t>(GDT_Float32);
	uint32_t gdalf64 = static_cast<uint32_t>(GDT_Float64);
	uint32_t gdali8 = static_cast<uint32_t>(GDT_Int8);
	uint32_t gdali16 = static_cast<uint32_t>(GDT_Int16);
	uint32_t gdali32 = static_cast<uint32_t>(GDT_Int32);

	//create band breaks vectors for float if required
	std::vector<std::vector<float>> fBandBreaks;
	if (dataPixelType == GDT_Float32) {
		fBandBreaks.resize(bandBreaks.size());
		for (size_t b = 0; i < bandBreaks.size(); b++) {
			fBandBreaks[b].resize(bandBreaks[b].size());
			for (size_t i = 0; i < bandBreaks[b].size(); i++) {
				fBandBreaks[b][i] = static_cast<float>(bandBreaks[b][i]);
			}
		}
	}

	//call the correct templated function depending on pixel types
	if (map) {
		switch (dataBandType | (stratBandType << 16)) {
			case (gdalf32 | (gdali8 << 16)):
				stratify<float, int8_t>(dataBands, stratBands, dataPixelType, stratPixelType,
							height, width, fBandBreaks, multipliers, threads);
				break;
			case (gdalf32 | (gdali16 << 16)):
				stratify<float, int16_t>(dataBands, stratBands, dataPixelType, stratPixelType,
							height, width, fBandBreaks, multipliers, threads);
				break;
			case (gdalf32 | (gdali32 << 16)):
				stratify<float, int32_t>(dataBands, stratBands, dataPixelType, stratPixelType,
							height, width, fBandBreaks, multipliers, threads);
				break;
			case (gdalf64 | (gdali8 << 16)):
				stratify<double, int8_t>(dataBands, stratBands, dataPixelType, stratPixelType,
							height, width, bandBreaks, multipliers, threads);
				break;
			case (gdalf64 | (gdali16 << 16)):
				stratify<double, int16_t>(dataBands, stratBands, dataPixelType, stratPixelType,
							height, width, bandBreaks, multipliers, threads);
				break;
			case (gdalf64 | (gdali32 << 16)):
				stratify<double, int32_t(dataBands, stratBands, dataPixelType, stratPixelType,
							height, width, bandBreaks, multipliers, threads);
				break;
			default:
				throw std::runtime_error("invalid combination of data band & strat band types. This is a bug.");
		}
	}
	else {	
			
	}

	//CLEAN UP allocated VECTOR!!!

	//close and add all of the VRT sub datasets as bands
	if (isVRTDataset) {
		for (size_t band = 0; band < VRTBandInfo.size(); band++) {
			GDALClose(VRTBandInfo[band].p_dataset);
			helper::addBandToVRTDataset(p_dataset, stratBands[band], VRTBandInfo[band]);
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
	       	new raster::GDALRasterWrapper(p_dataset) :
		new raster::GDALRasterWrapper(p_dataset, buffers);

}

} //namespace breaks
} //namespace sgs
