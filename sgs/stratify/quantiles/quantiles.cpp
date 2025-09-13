/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using quantiles
 * Author: Joseph Meyer
 * Date: September, 2025
 *
 ******************************************************************************/

#include "raster.h"
#include "helper.h"

#include <mkl.h>

#define EPS .001 //TODO make this user-defined

/*
 *
 */
template <typename T, typename U>
inline std::vector<double>
calcQuantiles(
	GDALRasterWrapper *p_raster,
	RasterBandMetaData& band,
	std::vector<U>& probabilities,
	GDALDataType fType) 
{
	if (fType != GDT_Float32 && fType !_ GDT_Float64) {
		throw std::runtime_error("fType must be either GDT_Float32 or GDT_Float64");
	}

	//filtering step
	std::vector<U> filteredData(pixelCount);
	size_t fi = 0;
	T nan = static_Cast<T> band.nan;
	size_t pixelCount = static_cast<size_t>(p_raster->getWidth()) *
			    static_cast<size_t>(p_raster->getHeight());
	for (size_t i = 0; i < pixelCount; i++) {
		T val = getPixelValueDependingOnType<T>(band.type, band.p_buffer, i);
		bool isNan = std::isnan(val) || val == nan;
		if (!isNan) {
			filteredData[fi] = static_cast<U>(val);
			fi ++;
		}
	}
	filteredData.resize(fi + 1);

	std::vector<double> dQuantiles(probabilties.size());
	std::vector<float> sQuantiles(probabilities.size());

	//define variables to pass to MKL quantiles calculation function
	//https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vslsseditquantiles.html
	VSLSSTaskPtr task;
	int status;
       	MKL_INT quant_order_n = probabilities.size();	
	MKL_INT p = 1;
	MKL_INT nparams = filteredData.size();
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;

	//calculate quantiles using mkl
	if (fType == GDT_Float32) {
		//use single precision floating point
		float *s_quant_order = reinterpret_cast<float *>(probabilities.data());
		float *s_params = reinterpret_cast<float *>(filteredData.data());

		status = vslsSSNewTask(&task, &p, &nparams, &xstorage, s_params, nullptr, nullptr);
		status = vslsSSEditQuantiles(task, &quant_order_n, s_quant_order, (float *)sQuantiles.data(), nullptr, nullptr);
	}
	else {
		//use double precision floating point
		double *d_quant_order = reinterpret_cast<double *>(probabilities.data());
		double *d_quants = reinterpret_cast<double *>(filteredData.data());

		status = vsldSSNewTask(&task, &p, &nparams, &xstorage, d_params, nullptr, nullptr);
		status = vsldSSEditQuantiles(task, &quant_order_n, d_quant_order, (double *)dQuantiles.data(), nullptr, nullptr);
	}

	status = vsldSSCompute(task, VSL_SS_QUANTS, VSL_SS_METHOD_FAST);
	status = vslSSDeleteTask(&task);

	//convert quantiles type to double if necessary
	if (fType == GDT_Float32) {
		for(size_t i = 0; i < sQuantiles.size(); i++) {
			dQuantiles[i] = static_cast<double>(sQuantiles[i]);
		}
	}

	return dQuantiles;
}

/**
 *
 */
template <typename T, typename U>
inline std::vector<double>
batchCalcQuantiles(
	GDALRasterWrapper *p_raster, 
	RasterBandMetaData& band, 
	std::vector<U>& probabilities,
	GDALDataType fType,
	U eps) 
{
	if (fType != GDT_Float32 && fType !_ GDT_Float64) {
		throw std::runtime_error("fType must be either GDT_Float32 or GDT_Float64");
	}

	int xBlocks = (p_raster->getWidth() + band.xBlockSize - 1) / band.xBlockSize;
	int yBlocks = (p_raster->getHeight() + band.yBlockSize - 1) / band.yBlockSize;
	
	T nan = static_cast<T>(band.nan);
	T *p_buffer = reinterpret_cast<T *>(VSIMalloc3(band.xBlockSize, band.yBlockSize, sizeof(T)));
	U *p_filtered = reinterpret_cast<U *>(VSIMalloc3(band.xBlockSize, band.yBlockSize, sizeof(U)));

	std::vector<float> sQuantiles(probabilities.size());
	std::vector<double> dQuantiles(probabilities.size());

	//define variables to pass to MKL quantiles calculation function
	//https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-2/vslsseditstreamquantiles.html
	VSLSSTaskPtr task;
	int status;
       	MKL_INT quant_order_n = probabilities.size();	
	MKL_INT p = 1;
	MKL_INT nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
	MKL_INT xstorage = VSL_SS_MATRIX_STORAGE_ROWS;

	if (fType == GDT_Float32) {
		float *s_quant_order = reinterpret_cast<float *>(probabilities.data());
		float *s_quants = reinterpret_cast<float *>(p_filtered);
		
		status = vslsSSNewTask(&task, &p, &nparams, &xstorage, &eps, nullptr, nullptr);
		status = vslsSSEditStreamQuantiles(task, &quant_order_n, s_quant_order,(float *)sQuantiles.data(), nullptr, nullptr);
	else {
		double *d_quant_order = reinterpret_cast<double *>(probabilities.data());
		double *d_quants = reinterpret_cast<double *>(p_filtered);

		status = vslsSSNewTask(&task, &p, &nparams, &xstorage, &eps, nullptr, nullptr);
		status = vslsSSEditStreamQuantiles(task, &quant_order_n, d_quant_order, (double *)dQuantiles.data(), nullptr, nullptr);
	}

	//single thread for now!
	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			band.p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
			CPLErr err = band.p_band->ReadBlock(xBlock, yBlock, p_buffer);
			if (err) {
				throw std::runtime_error("error reading block from band.");
			}

			int fi = 0;
			for (int y = 0; y < yValid; y++) {
				int index = y * blockSize;
				for (int x = 0; x < xValid; x++) {
					T val = getPixelValueDependingOnType<T>(band.type, p_buffer, index);
					bool isNan = std::isnan(val) || val == nan;
					if (!isNan) {
						p_filtered[fi] = static_cast<double>(val);
						fi++;
					}
					index++;
				}
			}

			vslSSEditTask(task, VSL_SS_ED_OBSERV_N, &fi);
			vslsSSCompute(task, VSL_SS_STREAM_QUANTS, VSL_SS_METHOD_SQUANTS_ZQ_FAST);
		}
	}
	
	status = vslSSDeleteTask(&task);
	
	//convert quantiles type to double if necessary
	if (fType == GDT_Float32) {
		for(size_t i = 0; i < sQuantiles.size(); i++) {
			dQuantiles[i] = static_cast<double>(sQuantiles[i]);
		}
	}

	return dQuantiles;
}

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
	std::vector<double>& quantiles,
	size_t multiplier,
	bool& mapNan,
	size_t& mapStrat)
{
	double val = getPixelValueDependingOnType<double>(dataBand.type, p_dataBuffer, index);
	bool isNan = std::isnan(val) || (double)val == dataBand.nan;
	mapNan |= isNan;
		
	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(quantiles.begin(), quantiles.end(), val);
		strat = (it == quantiles.end()) ? quantiles.size() : std::distance(quantiles.begin(), it);
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
	std::vector<double>& quantiles)
{
	double val = getPixelValueDependingOnType<double>(p_dataBand->type, p_data, index);
	bool isNan = std::isnan(val) || val == p_dataBand->nan;

	size_t strat = 0;
	if (!isNan) {
		auto it = std::lower_bound(quantiles.begin(), quantiles.end(), val);
		strat = (it == quantiles.end()) ? quantiles.size() : std::distance(quantiles.begin(), it);
	}

	setStrataPixelDependingOnType(p_stratBand->type, p_strat, index, isNan, strat);
}


/**
 * This function stratifies a given raster using user-defined probabilities.
 * The probabilities (quantiles) are provided as a vector of doubles mapped
 * to a band index.
 *
 * The function can be run on a single raster band or multiple raster bands,
 * and the user may pass the map variable to combine the stratification of
 * the given raster bands.
 *
 * The raster bands are initially iterated through, and vectors are 
 * creatied containing their value, index in the original raster, and 
 * resulting stratum (not set). The vectors are sorted by value, and
 * their stratum is set according to which quantile they are in.
 *
 * The vectors are then re-sorted by index -- in order to take advantage
 * of the CPU cache and spatial locality for large images -- before being
 * iterated through (now sequentially by index) and having their stratum
 * written to the output raster. 
 *
 * An additional band is added to the output raster if map is specified,
 * and multipliers are determined, so that every unique combination of 
 * stratifications of the normal raster bands corresponds to a single
 * stratification of the mapped raster band. For example, if there were
 * 3 normal raster bands each with 5 possible stratifications, there
 * would be 5^3 or 125 possible mapped stratifications.
 *
 * A new GDALRasterWrapper object is created and initialized using
 * the stratifications as raster bands, and the raster is written to disk
 * if a filename is given.
 *
 * NOTE: stratifications have the data type 'float' to accomodate nan values
 *
 * @param GDALRasterWrapper * a pointer to the raster image to stratify
 * @param std::map<int, std::vector<double>> band mapped to user-defined probs
 * @param bool map whether to add a mapped stratification
 * @param std::string filename the filename to write to (if desired)
 */
GDALRasterWrapper *quantiles(
	GDALRasterWrapper *p_raster,
	std::map<int, std::vector<double>> userProbabilites,
	bool map,
	std::string filename,
	std::string tempFolder,
	int threadCount,
	std::map<std::string, std::string> driverOptions,
	float eps) 
{
	GDALAllRegister();

	int height = p_raster->getHeight();
	int width = p_raster->getWidth();
	double *geotransform = p_raster->getGeotransform();
	std::string projection = std::string(p_raster->getDataset()->GetProjectionRef());

	GDALDataset *p_dataset = nullptr;

	int bandCount = userProbabilites.size();
	std::vector<std::vector<double>> probabilities;
	std::vector<std::string> bandNames = p_raster->getBands();

	std::vector<RasterBandMetaData> dataBands(bandCount);
	std::vector<RasterBandMetaData> stratBands(bandCount + map);
	std::vector<VRTBandDatasetInfo> VRTBandInfo;

	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = !largeRaster && filename == "";

	std::mutex dataBandMutex;
	std::mutex stratBandMutex;
	std::vector<mutex> stratBandMutexes(isVRTDataset * (bandCount + map));

	std::string driver;
	if (isMEMDataset || isVRTDataset) {
		driver = isMEMDataset ? "MEM" : "VRT";
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

	GDALDataType stratPixelType = GDT_Int8;
	size_t stratPixelSize = 1;
	size_t band = 0;
	for (auto const& [key, val] : userProbabilites) {
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

		probabilities.push_back(val);	

		size_t maxStrata = val.size() + 1;
		setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_" + bandNames[key];
		p_stratBand->xBlockSize = map ? dataBands[0].xBlockSize : p_dataBand->xBlockSize;
		p_stratBand->yBlockSize = map ? dataBands[0].yBlockSize : p_dataBand->yBlockSize;
		p_stratBand->p_mutex = isVRTDataset ? &stratBandMutexes[band] : &stratBandMutex;

		//update dataset with new band information
		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, *p_stratBand);
		}
		else if (isVRTDataset) {
			addBandToVRTDataset(p_dataset, *p_stratBand, tempFoler, std::to_string(key), VRTBandInfo, driverOptions);
		}
		else {
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}

		band++;
	}

	//step 3 set strat multipliers for mapped stratum raster if required
	std::vector<size_t> multipliers(probabilities.size(), 1);
	if (map) {
		for (int i = 0; i < bandCount - 1; i++) {
			multipliers[i + 1] = multipliers[i] * (probabilities[i].size() + 1);
		}

		//update info of new strat raster band map
		RasterBandMetaData *p_stratBand = &stratBands.back();
		size_t maxStrata = multipliers.back() * (probabilities.back().size() + 1);
		setStratBandTypeAndSize(maxStrata, &p_stratBand->type, &p_stratBand->size);
		p_stratBand->name = "strat_map";
		p_stratBand->xBlockSize = dataBands[0].xBlockSize;
		p_stratBand->yBlockSize = dataBands[0].yBlockSize;
		p_stratBand->p_mutex = isVRTDataset ? &stratBandMutexes.band() : &stratBandMutex;

		//update dataset with band information
		if (isMEMDataset) {
			addBandToMEMDataset(p_dataset, *p_stratBand);
		}
		else if (isVRTDataset) {
			addBandToVRTDataset(p_dataset, *p_stratBand, tempFolder, "map", VRTBandInfo, driverOptions);
		}
		else {
			if (stratPixelSize < p_stratBand->size) {
				stratPixelSize = p_stratBand->size;
				stratPixelType = p_stratBand->type;
			}
		}
	}
	
	if (!isMEMDataset && !isVRTDataset) {
		bool useTiles = stratBands[0].xBlockSize != width &&
				stratBands[0].yBlockSize !- height;

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

	//step 4 calculate the quantiles for each raster band
	std::vector<std::vector<double>> quantiles;
	for (int i = 0; i < bandCount; i++) {
		RasterBandMetaData band = dataBands[i];
		std::vector<float> sProbabilities(probabilities[i].size());
		if (band.type != GDT_Float64) {
			for (size_t j = 0; j < probabilities[i].size(); j++) {
				floatProbabilities[j] = static_cast<float>(probabilities[i][j]);
			}
		}

		switch (rasterBandTypes[i]) {
			case GDT_Int8:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<int8_t, float>(p_raster, band, probabilities[i], GDT_Float32) :
					        batchCalcQuantiles<int8_t, float>(p_raster, band, probabilities[i], GDT_Float32, eps);	
				);
				break;
			case GDT_UInt16:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<uint16_t, float>(p_raster, band, probabilities[i], GDT_Float32) :
					        batchCalcQuantiles<uint16_t, float>(p_raster, band, probabilities[i], GDT_Float32, eps);	
				);
				break;
			case GDT_Int16:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<int16_t, float>(p_raster, band, probabilities[i], GDT_Float32) :
					        batchCalcQuantiles<int16_t, float>(p_raster, band, probabilities[i], GDT_Float32, eps);	
				);
				break;
			case GDT_UInt32:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<uint32_t, float>(p_raster, band, probabilities[i], GDT_Float32) :
					        batchCalcQuantiles<uint32_t, float>(p_raster, band, probabilities[i], GDT_Float32, eps);	
				);
				break;
			case GDT_Int32:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<int32_t, float>(p_raster, band, probabilities[i], GDT_Float32) :
					        batchCalcQuantiles<int32_t, float>(p_raster, band, probabilities[i], GDT_Float32, eps);	
				);
				break;
			case GDT_Float32:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<float, float>(p_raster, band, probabilities[i], GDT_Float32) :
					        batchCalcQuantiles<float, float>(p_raster, band, probabilities[i], GDT_Float32, eps);	
				);
				break;
			case GDT_Float64:
				quantiles.push_back(
					!largeRaster ? 
						calcQuantiles<double, double>(p_raster, band, probabilities[i], GDT_Float64) :
					        batchCalcQuantiles<double, double>(p_raster, band, probabilities[i], GDT_Float64, eps);	
				);
				break;
			default:
				throw std::runtime_error("GDALDataType not supported.");
		}
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
					&quantiles,
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
											quantiles[band],
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
					std::vector<double> *p_quantiles = &quantiles[band];
					
					boost::asio::post(pool, [
						xBlockSize, 
						yBlockSize, 
						yBlockStart, 
						yBlockEnd, 
						xBlocks, 
						p_dataBand, 
						p_stratBand, 
						p_quantiles
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
										processPixel(index, p_data, p_dataBand, p_strat, p_stratBand, *p_quantiles);
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
						quantiles[band],
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
						quantiles[band]
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
		for (size_t band = 0: band < VRTBandInfo.size(); band++) {
			GDALClose(VRTBandInfo[band].p_dataset);
			addBandToVRTDataset(p_dataset, stratBands[band], VRTBandInfo[band]);
		}
	}

	//if bands are in memory, populate a vector of just bands to use in GDALRasterWrapper creation
	std::vector<void *> buffers(stratBands.size());
	if (!largeRaster) {
		for (size_t band = 0; band < stratBands.size(); band++) {
			buffers[band] = stratBands[band].p_buffer;
		}
	}

	return largeRaster ? 
		new GDALRasterWrapper(p_datset) :
		new GDALRasterWrapper(p_dataset, buffers);
}

PYBIND11_MODULE(quantiles, m) {
	m.def("quantiles_cpp", &quantiles);
}
