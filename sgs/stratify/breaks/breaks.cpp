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
	size_t bandCount,
	std::vector<GDALDataType>& rasterBandTypes,
	std::vector<void *>& rasterBandBuffers,
	std::vector<GDALDataType>& stratBandTypes,
	std::vector<void *>& stratBandBuffers,
	std::vector<size_t>& bandStratMultipliers,
	std::vector<std::vector<double>>& bandBreaks,
	std::vector<double> noDataValues)
{
	size_t mappedStrat = 0;
	bool mapNan = false;
	for (size_t i = 0; i < bandCount; i++) {
		double val = getPixelValueDependingOnType<double>(
			rasterBandTypes[i],
			rasterBandBuffers[i],
			index
		);
		bool isNan = std::isnan(val) || (double)val == noDataValues[i];
		mapNan |= isNan;
		
		//calculate strata value if not nan
		size_t strat = 0;
		if (!isNan) {
			std::vector<double> curBandBreaks = bandBreaks[i];
			auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
			strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
		}

		setStrataPixelDependingOnType(
			stratBandTypes[i],
			stratBandBuffers[i],
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
		stratBandTypes.back(),
		stratBandBuffers.back(),
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
	GDALDataType rasterBandType,
	void *rasterBandBuffer,
	GDALDataType stratBandType,
	void *stratBandBuffer,
	std::vector<double>& bandBreaks,
	double noDataValue)
{
	double val = getPixelValueDependingOnType<double>(
		rasterBandType,
		rasterBandBuffer,
		index
	);
	bool isNan = std::isnan(val) || (double)val == noDataValue;
		

	//calculate strata value if not nan
	size_t strat = 0;
	if (!isNan) {
		std::vector<double> curBandBreaks = bandBreaks;
		auto it = std::lower_bound(curBandBreaks.begin(), curBandBreaks.end(), val);
		strat = (it == curBandBreaks.end()) ? curBandBreaks.size() : std::distance(curBandBreaks.begin(), it);
	}

	setStrataPixelDependingOnType(
		stratBandType,
		stratBandBuffer,
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
	double *geotransform = p_raster->getGeotransform;
	std::string projection = std::string(p_raster->getDataset()->GetProjectionRef());

	GDALDataset *p_dataset = nullptr;

	size_t bandCount = breaks.size();
	std::vector<std::vector<double>> bandBreaks;

	std::vector<GDALRasterBand *> rasterBands;
	std::vector<void *> rasterBandBuffers;
	std::vector<GDALDataType> rasterBandTypes;
	std::vector<size_t> rasterBandTypeSizes;

	std::vector<GDALRasterBand *> stratBands;
	std::vector<void *> stratBandBuffers;
	std::vector<GDALDataType> stratBandTypes;
	std::vector<size_t> stratBandTypeSizes;

	std::vector<std::string> bandNames = p_raster->getBands();
	std::vector<double> noDataValues;

	bool isMEMDataset = !largeRaster && filename != "";
	bool isVRTDataset = largeRaster && filename != "";

	std::string driver;
	if (isMemDataset || isVRTDataset) {
		std::string driver = isMemDataset ? "MEM" : "VRT";
		p_dataset = createDataset(filename, driver, width, height, 0, GDT_Unknown, geotransform, projection);
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
	int firstXBlockSize = -1;
	int firstYBlockSize = -1;

	std::vector<GDALDataset *> VRTSubDatasets;

	GDALDataType largestStratType = GDT_Int8;
	size_t largestStratTypeSize = 1;
	for (auto const& [key, val] : breaks) {
		//update info of raster to stratify
		GDALRasterBand *p_band = p_raster->getRasterBand(key);
		rasterBands.push_back(p_band);
		rasterBandTypes.push_back(p_raster->getRasterBandType(key));
		rasterBandTypeSizes.push_back(p_raster->getRasterBandTypeSize(key));
		noDataValues.push_back(p_band->GetNoDataValue());
		if (!largeRaster) {
			rasterBandBuffers.push_back(p_raster->getRasterBandBuffer(key));
		}
	
		//determine block size to set the new strat raster as	
		int curXBlockSize, curYBlockSize;
		if (largeRaster) {
			rasterBands.back()->GetBlockSize(&curXBlockSize, &curYBlockSize);
			if (firstXBlockSize == -1) {
				firstXBlockSize = curXBlockSize;
			}
			if (firstYBlockSize == -1) {
				firstYBlockSize = curYBlockSize;
			}
		}

		//sort and add band breaks vector
		std::vector<double> valCopy = val; //have to create copy to alter the band breaks in iteration loop
		std::sort(valCopy.begin(), valCopy.end());
		bandBreaks.push_back(valCopy);

		//update info of new strat raster
		size_t pixelTypeSize = setStratBandType(val.size() + 1, stratBandTypes);
		stratBandTypeSizes.push_back(pixelTypeSize);
		
		if (isMEMDataset)
			addBandToMemDataset(
				stratBands,
				stratBandBuffers,
				p_dataset,
				stratBandTypes.back(),
				pixelTypeSize,
				"strat_" + bandNames[key]
			);
		}
		else if (isVRTDataset) {
			std::filesystem::path tmpPath = tempFolder;
			std::filesystem::path tmpName = "strat_breaks_" + bandNames[key]; //TODO potentially add timestamp to avoid same filenames
			tmpPath.append(tmpName);
			
			addBandToVRTDataset(
				stratBands,
				VRTSubDatasets,
				p_dataset,
				stratBandTypes.back(),
				"strat_" + bandNames[key]
				map ? firstXBlockSize : curXBlockSize,
				map ? firstYBlockSize : curYBlockSize,
				tmpPath.string(),
				geotransform,
				projection
			);			
		}
		else { //dataset driver is the one which corresponds to the filename extension
			if (largestStratTypeSize < pixelTypeSize) {
				largestStratTypeSize = pixelTypeSize;
				largestStratType = stratBandTypes.back();
			}	
		}
	}

	//step 2: set bandStratMultipliers and check max size if mapped stratification	
	std::vector<size_t> bandStratMultipliers(breaks.size(), 1);
	if (map) {
		//determine the stratification band index multipliers of the mapped band and error check maxes
		for (size_t i = 0; i < bandCount - 1; i++) {
			bandStratMultipliers[i + 1] = bandStratMultipliers[i] * (bandBreaks[i].size() + 1);
		}

		//update info of new strat raster map band
		size_t maxStrata = bandStratMultipliers.back() * (bandBreaks.back().size() + 1);
		size_t pixelTypeSize = setStratBandType(maxStrata, stratBandTypes);
		stratBandTypeSizes.push_back(pixelTypeSize);

		if (isMEMDataset)
			addBandToMemDataset(
				stratBands,
				stratBandBuffers,
				p_dataset,
				stratBandTypes.back(),
				pixelTypeSize,
				"strat_map"
			);
		}
		else if (isVRTDataset) {
			std::filesystem::path tmpPath = tempFolder;
			std::filesystem::path tmpName = "strat_breaks_map"; //TODO potentially add timestamp to avoid same filenames
			tmpPath.append(tmpName);

			addBandToVRTDataset(
				stratBands,
				VRTSubDatasets,
				p_dataset,
				stratBandTypes.back(),
				"strat_map",
				map ? firstXBlockSize : curXBlockSize,
				map ? firstYBlockSize : curYBlockSize,
				tmpPath.string(),
				geotransform,
				projection
			);			
		}
		else { //dataset driver is the one which corresponds to the filename extension
			if (largestStratTypeSize < pixelTypeSize) {
				largestStratTypeSize = pixelTypeSize;
				largestStratType = stratBandTypes.back();
			}	
		}
	}	

	//now we can create the dataset if the type is not MEM or VRT because we know the types of each band
	if (!isMEMDataset && !isVRTDataset) {
		p_dataset = createDataset(
			filename,
			driver,
			width,
			height,
			bandCount + static_cast<size_t>(map),
			largestStratType,
			geotransform,
			projection
		);

		//for each band, update bands, size, type, and potentially allocate a buffer
		for (size_t band = 0; band < bandCount + static_cast<size_t>(map); band++) {
			stratBandTypes[band] = largestStratType;
			stratBandTypesSizes[band] = largestStratTypeSize;
			stratBands.push_back(p_dataset->GetRasterBand(band + 1));
			if (!largeRaster) {
				void *p_data = VSIMalloc3(height, width, largestStraTypeSize);
				stratBandBuffers.push_back(p_data);
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
			int xBlockSize, yBlockSize, checkXBlockSize, checkYBlockSize;
			rasterBands[0]->GetBlockSize(&xBlockSize, &yBlockSize);

			for (size_t band = 0; band < bandCount; band++) {
				rasterBands[band]->GetBlockSize(&checkXBlockSize, &checkYBlockSize);
				useBlocksRaster.push_back(checkXBlockSize == xBlockSize && checkYBlockSize == yBlockSize);
				rasterBandBuffers.push_back(VSIMalloc3(xBlockSize, yBlockSize, rasterBandTypeSizes[band]));
				
				stratBands[band]->GetBlockSize(&checkXBlockSize, &checkYBlockSize);
				useBlocksStrat.push_back(checkXBlockSize == xBlockSize && checkYBlockSize == yBlockSize);
				stratBandBuffers.push_back(VSIMalloc3(xBlockSize, yBlockSize, stratBandTypeSizes[band]));
			}	

			//mapped strat raster block size calculation
			stratBands[bandCount]->GetBlockSize(&checkXBlockSize, &checkYBlockSize);
			useBlocksStrat.push_back(checkXBlockSize == xBlockSize && checkYBlockSize == yBlockSize);
			stratBandBuffers.push_back(VSIMalloc3(xBlockSize, yBlockSize, stratBandTypeSizes.back()));

			int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
			int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;

			for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
				for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
					int xValid, yValid;
					rasterBands[0]->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);

					//read raster band data into band buffers
					for (size_t band = 0; band < bandCount; band++) {
						if (useBlocksRaster[band]) {
							rasterBands[band]->ReadBlock(
								xBlock, 
								yBlock, 
								rasterBandBuffers[band]
							);
						}
						else {
							rasterBands[band]->RasterIO(
								GF_Read,
								xBlock * xBlockSize,
								yBlock * yBlockSize,
								xValid,
								yValid,
								rasterBandBuffers[band],
								xBlockSize,
								yBlockSize,
								rasterBandTypes[band],
								0,
								(xBlockSize - xValid) * rasterBandTypeSizes[band]
							);
						}
					}

					//process blocked band data
					for (int y = 0; y < yValid; y++) {
						for (int x = 0; x < xValid; x++) {
							size_t index = static_cast<size_t>(x + y * xBlockSize);
							processMapPixel(
								index,
								bandCount,
								rasterBandTypes,
								rasterBandBuffers,
								stratBandTypes,
								stratBandBuffers,
								bandStratMultipliers,
								bandBreaks,
								noDataValues
							);
						}
					}
					
					//write strat band data
					for (size_t band = 0; band <= bandCount; band++) {
						if (useBlocksStrat[band]) {
							stratBands[band]->WriteBlock(
								xBlock, 
								yBlock, 
								stratBandBuffers[band]
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
					}
				}
			}	
		}
		else {
			std::cout << "HERE!" << std::endl;
			void *p_rast = nullptr;
			void *p_strat = nullptr;
			for (size_t band = 0; band < bandCount; band++) {
				int xBlockSize, yBlockSize, stratXBlockSize, stratYBlockSize;
				rasterBands[band]->GetBlockSize(&xBlockSize, &yBlockSize);
				stratBands[band]->GetBlockSize(&stratXBlockSize, &stratYBlockSize);
				
				std::cout << "xBlockSize: " << xBlockSize << std::endl;
				std::cout << "yBlockSize: " << yBlockSize << std::endl;

				if (xBlockSize != stratXBlockSize || yBlockSize != stratYBlockSize) {
					throw std::runtime_error("block size error -- should not be here!!!");
				}

				if (band == 0) {
					p_rast = VSIMalloc3(xBlockSize, yBlockSize, rasterBandTypeSizes[band]);
					p_strat = VSIMalloc3(xBlockSize, yBlockSize, stratBandTypeSizes[band]);
				}
				else {
					p_rast = VSIRealloc(p_rast, xBlockSize * yBlockSize * rasterBandTypeSizes[band]);
					p_strat = VSIRealloc(p_strat, xBlockSize * yBlockSize * stratBandTypeSizes[band]);
				}

				int xBlocks = (p_raster->getWidth() + xBlockSize - 1) / xBlockSize;
				int yBlocks = (p_raster->getHeight() + yBlockSize - 1) / yBlockSize;

				std::cout << "xBlocks: " << xBlocks << std::endl;
				std::cout << "yBlocks: " << yBlocks << std::endl;
				std::cout << "total blocks: " << (xBlocks * yBlocks) << std::endl;

				for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
					if (yBlock == 10) {
						throw std::runtime_error("debugging stop -- remove later");
					}
					for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
						rasterBands[band]->ReadBlock(xBlock, yBlock, p_rast);

						int xValid, yValid;
						rasterBands[band]->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);

						for (int y = 0; y < yValid; y++) {
							for (int x = 0; x < xValid; x++) {
								size_t index = static_cast<size_t>(x + y * xBlockSize);
								processPixel(
									index,
									rasterBandTypes[band],
									p_rast,
									stratBandTypes[band],
									p_strat,
									bandBreaks[band],
									noDataValues[band]
								);
							}
						}

						stratBands[band]->WriteBlock(xBlock, yBlock, p_strat);
					}
				}
			}

			VSIFree(p_rast);
			VSIFree(p_strat);
		}
	}
	else {
		size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) * static_cast<size_t>(p_raster->getHeight());
		if (map) {
			for (size_t i = 0; i < pixelCount; i++) {
				processMapPixel(
					i,
					bandCount,
					rasterBandTypes,
					rasterBandBuffers,
					stratBandTypes,
					stratBandBuffers,
					bandStratMultipliers,
					bandBreaks,
					noDataValues
				);		
			}
		}
		else {
			for (size_t band = 0; band < bandCount; band++) {
				for (size_t i = 0; i < pixelCount; i++) {
					processPixel(
						i, 
						rasterBandTypes[band],
						rasterBandBuffers[band],
						stratBandTypes[band],
						stratBandBuffers[band],
						bandBreaks[band],
						noDataValues[band]
					);
				}
			}
		}

		if (!isVRTDataset && !isMEMDataset) {
			for (size_t i = 0; i < stratBands.size(); i++) {
				GDALRasterBand *p_band = stratBands[i];
				void *p_data = stratBandBuffers[i];

				p_band->RasterIO(
					GF_Write,
					0,
					0,
					width,
					height,
					p_data,
					width,
					height, 
					largestStratType,
					0,
					0
				);
			}
		}
	}

	//free allocated band data
	if (map && largeRaster) {
		for (size_t i = 0; i < stratBandBuffers.size(); i++) {
			VSIFree(stratBandBuffers[i]);
		}
	}

	//step 4: create GDALRasterWrapper object from bands
	//this dynamically-allocated object will be cleaned up by python (TODO I hope...)
	GDALRasterWrapper *p_stratRaster = largeRaster ?
	       	new GDALRasterWrapper(p_dataset) :
		new GDALRasterWrapper(p_dataset, stratBandBuffers);

	//step 5: close all of the VRT sub datasets
	if (isVRTDataset) {
		for (size_t i = 0; i < VRTSubDatasets.size(); i++) {
			GDALClose(VRTSubDatasets[i]);
		}
	}
	
	return p_stratRaster;
}

PYBIND11_MODULE(breaks, m) {
	m.def("breaks_cpp", &breaks);
}
