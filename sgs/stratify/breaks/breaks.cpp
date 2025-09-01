/******************************************************************************
 *
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using breaks
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#include <iostream>

#include "raster.h"
#include "helper.h"

/**
 *
 */
inline void processMapPixel(
	size_t index,
	std::vector<RasterBandMetaData>& dataBands,
	std::vector<RasterBandMetaData>& stratBands,
	std::vector<size_t>& bandStratMultipliers,
	std::vector<std::vector<double>>& bandBreaks)
{
	size_t mappedStrat = 0;
	bool mapNan = false;
	for (size_t i = 0; i < dataBands.size(); i++) {
		double val = getPixelValueDependingOnType<double>(
			dataBands[i].type,
			dataBands[i].p_buffer,
			index
		);
		bool isNan = std::isnan(val) || (double)val == dataBands[i].nan;
		mapNan |= isNan;
		
		//calculate strata value if not nan
		size_t strat = 0;
		if (!isNan) {
			std::vector<double> curBandBreaks = bandBreaks[i];
			auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
			strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
		}

		setStrataPixelDependingOnType(
			stratBands[i].type,
			stratBands[i].p_buffer,
			index,
			isNan,
			strat
		);

		//adjust mappedStrat as required
		if (!mapNan) {
			mappedStrat += strat * bandStratMultipliers[i];
		}
	}
		
	//assign mapped value
	setStrataPixelDependingOnType(
		stratBands.back().type,
		stratBands.back().p_buffer,
		index,
		mapNan,
		mappedStrat
	);
}	

/**
 *
 */
inline void processPixel(
	size_t index,
	RasterBandMetaData dataBand,
	RasterBandMetaData stratBand,
	std::vector<double>& bandBreaks)
{
	double val = getPixelValueDependingOnType<double>(
		dataBand.type,
		dataBand.p_buffer,
		index
	);
	bool isNan = std::isnan(val) || (double)val == dataBand.nan;
		

	//calculate strata value if not nan
	size_t strat = 0;
	if (!isNan) {
		std::vector<double> curBandBreaks = bandBreaks;
		auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
		strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
	}

	setStrataPixelDependingOnType(
		stratBand.type,
		stratBand.p_buffer,
		index,
		isNan,
		strat
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
	std::vector<std::string> newBandNames;

	std::vector<RasterBandMetaData> dataBands;
	std::vector<RasterBandMetaData> stratBands;
	std::vector<VRTBandDatasetInfo> VRTBandInfo;

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";

	std::string driver;
	if (isMEMDataset || isVRTDataset) {
		std::string driver = isMEMDataset ? "MEM" : "VRT";
		p_dataset = createDataset(filename, driver, width, height, 0, {}, GDT_Unknown, geotransform, projection, nullptr);
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();

		if (extension != ".tif") {
			throw std::runtime_error("sgs only supports .tif files right now");
		}

		driver = "Gtiff";
	}

	//step 1: allocate, read, and initialize raster data and breaks information

	GDALDataType stratPixelType = GDT_Int8;
	size_t stratPixelSize = 1;
	for (auto const& [key, val] : breaks) {
		//update info of raster to stratify
		GDALRasterBand *p_band = p_raster->getRasterBand(key);
		
		RasterBandMetaData dataBand;
		dataBand.p_band = p_band;
		dataBand.type = p_raster->getRasterBandType(key);
		dataBand.size = p_raster->getRasterBandTypeSize(key);
		dataBand.p_buffer = largeRaster ? nullptr : p_raster->getRasterBandBuffer(key);
		dataBand.nan = p_band->GetNoDataValue();
		p_band->GetBlockSize(&dataBand.xBlockSize, &dataBand.yBlockSize);
		dataBands.push_back(dataBand);

		//sort and add band breaks vector
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);

		//update info of new strat raster
		RasterBandMetaData stratBand;
		size_t maxStrata = val.size() + 1;
		setStratBandTypeAndSize(maxStrata, &stratBand.type, &stratBand.size);
		stratBand.name = "strat_" + bandNames[key];
		stratBand.xBlockSize = map ? dataBands[0].xBlockSize : dataBand.xBlockSize;
		stratBand.yBlockSize = map ? dataBands[0].yBlockSize : dataBand.yBlockSize;
		
		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, stratBand);
		}
		else if (isVRTDataset) {
			std::filesystem::path tmpPath = tempFolder;
			std::filesystem::path tmpName = "strat_breaks_" + std::to_string(key) + ".tif";
			tmpPath = tmpPath / tmpName;
		
			VRTBandDatasetInfo info;		
			info.filename = tmpPath.string();
			createVRTSubDataset(p_dataset, stratBand, info);
			VRTBandInfo.push_back(info);
		}
		else { //dataset driver is the one which corresponds to the filename extension
			newBandNames.push_back(stratBand.name);
			if (stratPixelSize < stratBand.size) {
				stratPixelSize = stratBand.size;
				stratPixelType = stratBand.type;
			}
		}

		stratBands.push_back(stratBand);
	}

	//step 2: set bandStratMultipliers and check max size if mapped stratification	
	std::vector<size_t> bandStratMultipliers(breaks.size(), 1);
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (size_t i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (bandBreaks[i].size() + 1);
		}

		//update info of new strat raster map band
		RasterBandMetaData stratBand;
		size_t maxStrata = bandStratMultipliers.back() * (bandBreaks.back().size() + 1);
		setStratBandTypeAndSize(maxStrata, &stratBand.type, &stratBand.size);
		stratBand.name = "strat_map";
		stratBand.xBlockSize = dataBands[0].xBlockSize;
		stratBand.yBlockSize = dataBands[0].yBlockSize;

		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, stratBand);
		}
		else if (isVRTDataset) {
			std::filesystem::path tmpPath = tempFolder;
			std::filesystem::path tmpName = "strat_breaks_map.tif";
			tmpPath = tmpPath / tmpName;

			VRTBandDatasetInfo info;
			info.filename = tmpPath.string();
			createVRTSubDataset(p_dataset, stratBand, info);
			VRTBandInfo.push_back(info);
		}
		else { //dataset driver is the one which corresponds to the filename extension
			newBandNames.push_back(stratBand.name);
			if (stratPixelSize < stratBand.size) {
				stratPixelType = stratBand.type;
				stratPixelSize = stratBand.size;
			}
		}

		stratBands.push_back(stratBand);
	}	

	//now we can create the dataset if the type is not MEM or VRT because we know the types of each band
	if (!isMEMDataset && !isVRTDataset) {
		char **papszOptions = nullptr;
		if (largeRaster) {
			papszOptions = CSLSetNameValue(
				papszOptions,
				"BLOCKXSIZE",
				std::to_string(stratBands[0].xBlockSize).c_str()
			);
			papszOptions = CSLSetNameValue(
				papszOptions,
				"BLOCKYSIZE",
				std::to_string(stratBands[0].xBlockSize).c_str()
			);
		}
		
		p_dataset = createDataset(
			filename,
			driver,
			width,
			height,
			bandCount + static_cast<size_t>(map),
			newBandNames,
			stratPixelType,
			geotransform,
			projection,
			papszOptions
		);

		//for each band, update bands, size, type, and potentially allocate a buffer
		for (size_t band = 0; band < bandCount + static_cast<size_t>(map); band++) {
			stratBands[band].size = stratPixelSize;
			stratBands[band].type = stratPixelType;
			stratBands[band].p_band = p_dataset->GetRasterBand(band + 1);
			if (!largeRaster) {
				stratBands[band].p_buffer = VSIMalloc3(height, width, stratPixelSize);
			}
		}
	}

	//TODO: multithread and consider cache thrashing
	//step 3: iterate through indices and update the stratified raster bands
	if (largeRaster) {
		if (map) {
			std::vector<bool> useBlocksRaster;
			std::vector<bool> useBlocksStrat;

			//use the first raster band to determine block size
			int xBlockSize = dataBands[0].xBlockSize;
			int yBlockSize = dataBands[0].yBlockSize;

			for (size_t band = 0; band < bandCount; band++) {
				useBlocksRaster.push_back(
					dataBands[band].xBlockSize == xBlockSize && 
					dataBands[band].yBlockSize == yBlockSize
				);
				dataBands[band].p_buffer = VSIMalloc3(xBlockSize, yBlockSize, dataBands[band].size);
				
				useBlocksStrat.push_back(
					stratBands[band].xBlockSize == xBlockSize && 
					stratBands[band].yBlockSize == yBlockSize
				);
				stratBands[band].p_buffer = VSIMalloc3(xBlockSize, yBlockSize, stratBands[band].size);
			}	

			//mapped strat raster block size calculation
			useBlocksStrat.push_back(
				stratBands.back().xBlockSize == xBlockSize && 
				stratBands.back().yBlockSize == yBlockSize
			);
			stratBands.back().p_buffer = VSIMalloc3(xBlockSize, yBlockSize, stratBands.back().size);

			int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
			int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;

			for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
				for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
					int xValid, yValid;
					dataBands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);

					//read raster band data into band buffers
					CPLErr err;
					for (size_t band = 0; band < bandCount; band++) {
						if (useBlocksRaster[band]) {
							err = dataBands[band].p_band->ReadBlock(
								xBlock, 
								yBlock, 
								dataBands[band].p_buffer
							);
						}
						else {
							err = dataBands[band].p_band->RasterIO(
								GF_Read,
								xBlock * xBlockSize,
								yBlock * yBlockSize,
								xValid,
								yValid,
								dataBands[band].p_buffer,
								xBlockSize,
								yBlockSize,
								dataBands[band].type,
								0,
								(xBlockSize - xValid) * dataBands[band].size
							);
						}
						if (err) {
							throw std::runtime_error("unable to read block from raster");
						}
					}

					//process blocked band data
					for (int y = 0; y < yValid; y++) {
						for (int x = 0; x < xValid; x++) {
							size_t index = static_cast<size_t>(x + y * xBlockSize);
							processMapPixel(index, dataBands, stratBands, bandStratMultipliers, bandBreaks);
						}
					}
					
					//write strat band data
					for (size_t band = 0; band <= bandCount; band++) {
						if (useBlocksStrat[band]) {
							err = stratBands[band].p_band->WriteBlock(
								xBlock, 
								yBlock, 
								stratBands[band].p_buffer
							);
						}
						else {
							throw std::runtime_error("should not be here!!!");
							//stratBands[band]->RasterIO(
							//	GF_Write,
							//	xBlock * xBlockSize,
							//	yBlock * yBlockSize,
							//	xValid,
							//	yValid,
							//	stratBandBuffers[band],
							//	xBlockSize,
							//	yBlockSize,
							//	stratBandTypes[band],
							//	0,
							//	(xBlockSize - xValid) * stratBandTypeSizes[band]
							//);
						}
						if (err) {
							throw std::runtime_error("error writing block of data to band.");
						}
					}
				}
			}	
		}
		else {
			std::cout << "HERE!" << std::endl;
			void *p_data = nullptr;
			void *p_strat = nullptr;
			for (size_t band = 0; band < bandCount; band++) {
				int xBlockSize = dataBands[band].xBlockSize;
				int yBlockSize = dataBands[band].yBlockSize;
				std::cout << "xBlockSize: " << xBlockSize << std::endl;
				std::cout << "yBlockSize: " << yBlockSize << std::endl;
					
				if (xBlockSize != stratBands[band].xBlockSize || 
				    yBlockSize != stratBands[band].yBlockSize) {
					throw std::runtime_error("block size error -- should not be here!!!");
				}

				if (band == 0) {
					p_data = VSIMalloc3(xBlockSize, yBlockSize, dataBands[band].size);
					p_strat = VSIMalloc3(xBlockSize, yBlockSize, stratBands[band].size);
				}
				else {
					p_data = VSIRealloc(p_data, xBlockSize * yBlockSize * dataBands[band].size);
					p_strat = VSIRealloc(p_strat, xBlockSize * yBlockSize * stratBands[band].size);
				}
				dataBands[band].p_buffer = p_data;
				stratBands[band].p_buffer = p_strat;


				int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
				int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;

				std::cout << "xBlocks: " << xBlocks << std::endl;
				std::cout << "yBlocks: " << yBlocks << std::endl;
				std::cout << "total blocks: " << (xBlocks * yBlocks) << std::endl;

				for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
					if (yBlock == 10) {
						throw std::runtime_error("debugging stop -- remove later");
					}
					CPLErr err;
					for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
						err = dataBands[band].p_band->ReadBlock(xBlock, yBlock, dataBands[band].p_buffer);
						if (err) {
							throw std::runtime_error("error reading block from band.");
						}

						int xValid, yValid;
						dataBands[band].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);

						for (int y = 0; y < yValid; y++) {
							for (int x = 0; x < xValid; x++) {
								size_t index = static_cast<size_t>(x + y * xBlockSize);
								processPixel(index, dataBands[band], stratBands[band], bandBreaks[band]);
							}
						}

						err = stratBands[band].p_band->WriteBlock(xBlock, yBlock, stratBands[band].p_buffer);
						if (err) {
							throw std::runtime_error("error writing block of data do band.");
						}
					}
				}
			}

			VSIFree(p_data);
			VSIFree(p_strat);
		}
	}
	else {
		size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) * static_cast<size_t>(p_raster->getHeight());
		if (map) {
			for (size_t i = 0; i < pixelCount; i++) {
				processMapPixel(i, dataBands, stratBands, bandStratMultipliers, bandBreaks);		
			}
		}
		else {
			for (size_t band = 0; band < bandCount; band++) {
				for (size_t i = 0; i < pixelCount; i++) {
					processPixel(i, dataBands[band], stratBands[band], bandBreaks[band]);
				}
			}
		}

		if (!isVRTDataset && !isMEMDataset) {
			for (size_t band = 0; band < stratBands.size(); band++) {
				stratBands[band].p_band->RasterIO(
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
