/******************************************************************************
 *
 * Project: sgs
 * Purpose: C++ implementation of PCA
 * Author: Joseph Meyer
 * Date: October, 2025
 *
 ******************************************************************************/

#include <iostream>

#include "helper.h"
#include "raster.h"

#include "oneapi/dal.hpp"

typedef oneapi::dal::homogen_table	DALHomogenTable;

/**
 * This struct contains the output eigenvectors and eigenvalues for
 * the principal components. It also contains the mean and standard
 * deviation for each raster band.
 *
 * The mean and standard deviation are used when writing the outputs
 * to both center and scale each band. The eigenvectors are then
 * used when calculating the output principal component values.
 */
template <typename T>
struct PCAResult {
	std::vector<std::vector<T>> eigenvectors;
	std::vector<T> eigenvalues;
	std::vector<double> means;
	std::vector<double> stdevs;
};

/**
 * This struct contains the intermediate values, as well as functions
 * for updating the intermediate values of the variance of a raster
 * band using Welfords method. The mean and standard deviation of
 * a raster band can be calculated afterwards without requiring the
 * whole raster to be in memory.
 *
 * Double precision is always used for higher precision of potentially small values.
 *
 * https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
 */
struct Variance {
	int64_t k;	//count
	double M = 0;	//running mean
	double S = 0;	//sum of squares
	double oldM = 0;
	
	inline void
	update(double x) {
		k++;
		oldM = M;

		//update running mean
		M = M + (x - M) / static_cast<double>(k);

		//update sum of squares
		S = S + (x - M) * (x - oldM);
	}

	inline double 
	getMean() {
		return M;
	}

	inline double 
	getStdev() {
		double variance = S / static_cast<double>(k);
		return std::sqrt(variance);
	}
};

/**
 * This function is used by the pca() function to calculate the principal component
 * eigenvectors and eigenvalues, along with the mean and standard deviation of each
 * input raster band. This function is used in the case where the input raster is
 * small, and can reasonably be expected to fit entirely into memory.
 *
 * First, the input raster bands are read into memory usign the GDALRasterBand 
 * RasterIO function. Bands are read into memory in a row-wise manor 
 * such that a row indicates a single pixel, and a column indicates a raster band.
 * This means that in between each pixel and the next, a gap must be left for the
 * remaining band values for that pixel index to be written to. This is done
 * using the nPixelSpace, and nLineSpace arguments of RasterIO.
 *
 * Second, each pixel is checked to ensure it isn't a nan pixel. Any pixel
 * containing a nan value in any band is overwritten completely with the
 * next not-nan pixel, the total number of not-nan pixels is stored as the
 * number of features.
 *
 * The mean, standard deviation are then calculated using Welfords method,
 * and the pca eigenvectors and eigenvalues are calculated using the oneDAL
 * library principal components functionality.
 *
 * A result containing the eigenvectors, eigenvalues, mean per band, and
 * standard deviation per band, is returned.
 */
template <typename T>
PCAResult<T>
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int width,
	int height,
	int nComp)
{	
	int bandCount = static_cast<int>(bands.size());
	T *p_data = reinterpret_cast<T *>(VSIMalloc3(width * height, bandCount, size));

	std::vector<Variance> bandVariances(bandCount);
	std::vector<T> noDataVals(bandCount);
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}

	//read full bands into p_data
	for (size_t i = 0; i < bands.size(); i++) {
		bands[i].p_band->RasterIO(
			GF_Read,
			0,
			0,
			width,
			height,
			(void *)((size_t)p_data + i * size),	
			width,
			height,
			type,
			size * bands.size(),
			size * bands.size() * width
		);	
	}

	//remove nodata values from p_data, ensuring data pixels are consecutive
	int nFeatures = 0;
	for (int i = 0; i < height * width; i++) {
		bool isNan = false;
		for (int b = 0; b < bandCount; b++) { 
			T val = p_data[i * bandCount + b];
			isNan = val == noDataVals[b] || std::isnan(val);
			if (isNan) {
				break;
			}
			p_data[nFeatures * bandCount + b] = val;
		}
		nFeatures += !isNan;
	}

	//update variance calculations
	for (int i = 0; i < nFeatures; i++) {
		for (int b = 0; b < bandCount; b++) {
			T val = p_data[i * bandCount + b];
			bandVariances[b].update(static_cast<double>(val));	
		}
	}

	//calculate pca
	DALHomogenTable table = DALHomogenTable::wrap<T>(p_data, nFeatures, bandCount, oneapi::dal::data_layout::row_major);
	const auto desc = oneapi::dal::pca::descriptor<float, oneapi::dal::pca::method::cov>().set_component_count(nComp).set_deterministic(true);
       	const auto result = oneapi::dal::train(desc, table);

	VSIFree(p_data);

	PCAResult<T> retval;
	auto eigenvectors = result.get_eigenvectors();
	auto eigenvalues = result.get_eigenvalues();
	int64_t eigRows = eigenvectors.get_row_count();
	int64_t eigCols = eigenvectors.get_column_count();

	oneapi::dal::row_accessor<const float> eigVecAcc {eigenvectors};
	auto eigVecBlock = eigVecAcc.pull({0, eigRows});

	oneapi::dal::row_accessor<const float> eigValAcc {eigenvalues};
	auto eigValBlock = eigValAcc.pull({0, 1});

	retval.eigenvectors.resize(eigRows);
	retval.eigenvalues.resize(eigRows);
	for (int64_t i = 0; i < eigRows; i++) {
		retval.eigenvalues[i] = static_cast<double>(eigValBlock[i]);

		retval.eigenvectors[i].resize(eigCols);
		for (int64_t j = 0; j < eigCols; j++) {
			retval.eigenvectors[i][j] = static_cast<double>(eigVecBlock[i * eigCols + j]);
		}
	}

	for (int b = 0; b < bandCount; b++) {
		retval.means.push_back(bandVariances[b].getMean());
		retval.stdevs.push_back(bandVariances[b].getStdev());
	}

	return retval;	
}

/**
 * This function is used by the pca() function to calculate the principal component
 * eigenvectors and eigenvalues, along with the mean and standard deviation of each
 * input raster band. This function is used in the case where the input raster is
 * large, will be processed in blocks.
 *
 * All of the blocks are iterated through, and within each iteration the following
 * is done:
 *
 * First, the input raster band blocks are read into memory using the GDALRasterBand 
 * RasterIO function. Bands are read into memory in a row-wise manor 
 * such that a row indicates a single pixel, and a column indicates a raster band.
 * This means that in between each pixel and the next, a gap must be left for the
 * remaining band values for that pixel index to be written to. This is done
 * using the nPixelSpace, and nLineSpace arguments of RasterIO.
 *
 * Second, each pixel is checked to ensure it isn't a nan pixel. Any pixel
 * containing a nan value in any band is overwritten completely with the
 * next not-nan pixel, the total number of not-nan pixels is stored as the
 * number of features.
 *
 * The mean, standard deviation are then updated using Welfords method,
 * and the pca eigenvectors and eigenvalues partial result are updated 
 * using the oneDAL library principal components functionality.
 *
 * once all blocks have been iterated through, the final resulting mean per band,
 * standard deviation per band, eigenvectors, and eigenvalues are calculated 
 * and returned.
 */
template <typename T>
PCAResult<T>
calculatePCA(
	std::vector<RasterBandMetaData>& bands,
	GDALDataType type,
	size_t size,
	int xBlockSize,
	int yBlockSize,
	int xBlocks,
	int yBlocks,
	int nComp)
{
	int bandCount = static_cast<int>(bands.size());
	T *p_data = reinterpret_cast<T *>(VSIMalloc3(xBlockSize * yBlockSize, bandCount, size));

	std::vector<Variance> bandVariances(bandCount);
	std::vector<T> noDataVals(bandCount);
	for (size_t i = 0; i < bands.size(); i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}

	const auto desc = oneapi::dal::pca::descriptor<float, oneapi::dal::pca::method::cov>().set_component_count(nComp).set_deterministic(true);
	oneapi::dal::pca::partial_train_result<> partial_result;

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			bands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
		
			//read bands into memory	
			for (size_t i = 0; i < bandCount; i++) {
				bands[i].p_band->RasterIO(
					GF_Read,
					xBlock * xBlockSize,
					yBlock * yBlockSize,
					xValid,
					yValid,
					(void *)((size_t)p_data + i * size),
					xValid,
					yValid,
					type,
					size * bands.size(),
					size * bands.size() * static_cast<size_t>(xBlockSize)
				);				
			}

			//remove nodata values
			int nFeatures = 0;
			for (int x = 0; x < xValid; x++) {
				for (int y = 0; y < yValid; y++) {
					bool isNan = false;
					for (int b = 0; b < bandCount; b++) {
						T val = p_data[((y * xBlockSize) + x) * bandCount + b];
						isNan = std::isnan(val) || val == noDataVals[b];
						if (isNan) {
							break;
						}
						p_data[nFeatures * bandCount + b] = val;
					}
					nFeatures += !isNan;
				}
			}

			//update variance calculations
			for (int i = 0; i < nFeatures; i++) {
				for (int b = 0; b < bandCount; b++) {
					T val = p_data[i * bandCount + b];
					bandVariances[b].update(static_cast<double>(val));	
				}
			}

			//calculate partial result
			DALHomogenTable table = DALHomogenTable(p_data, nFeatures, bandCount, [](const T*){}, oneapi::dal::data_layout::row_major);
			partial_result = oneapi::dal::partial_train(desc, partial_result, table);
		}
		
		if (yBlock % 1000 == 0) {
			std::cout << "completed y block " << yBlock << "/" << yBlocks << std::endl;
		}
	}

	auto result = oneapi::dal::finalize_train(desc, partial_result);
	
	VSIFree(p_data);

	PCAResult<T> retval;
	auto eigenvectors = result.get_eigenvectors();
	auto eigenvalues = result.get_eigenvalues();
	int64_t eigRows = eigenvectors.get_row_count();
	int64_t eigCols = eigenvectors.get_column_count();

	oneapi::dal::row_accessor<const float> eigVecAcc {eigenvectors};
	auto eigVecBlock = eigVecAcc.pull({0, eigRows});

	oneapi::dal::row_accessor<const float> eigValAcc {eigenvalues};
	auto eigValBlock = eigValAcc.pull({0, 1});

	retval.eigenvectors.resize(eigRows);
	retval.eigenvalues.resize(eigRows);
	for (int64_t i = 0; i < eigRows; i++) {
		retval.eigenvalues[i] = static_cast<T>(eigValBlock[i]);

		retval.eigenvectors[i].resize(eigCols);
		for (int64_t j = 0; j < eigCols; j++) {
			retval.eigenvectors[i][j] = static_cast<T>(eigValBlock[i * eigCols + j]);
		}
	}
	
	for (int b = 0; b < bandCount; b++) {
		retval.means.push_back(bandVariances[b].getMean());
		retval.stdevs.push_back(bandVariances[b].getStdev());
	}

	return retval;	
}


/**
 * This function is used to write the output principal components to a
 * raster dataset, after the eigenvectors and eigenvalues have already
 * been calculated for the input raster. This function is used in the
 * case where the raster is small, and would not be expected to cause
 * errors for being entirely in memory.
 *
 * First, the input raster bands are read into memory using the GDALRasterBand 
 * RasterIO function. Bands are read into memory in a row-wise manor 
 * such that a row indicates a single pixel, and a column indicates a raster band.
 * This means that in between each pixel and the next, a gap must be left for the
 * remaining band values for that pixel index to be written to. This is done
 * using the nPixelSpace, and nLineSpace arguments of RasterIO.
 *
 * Second, the pixels are iterated over. If any value in any band is a no-data
 * value, then nan is written. If none of the values for any band are nodata,
 * then the output pca value is calculated. A different function is called
 * depending on whether the data type is float (single precision), or double 
 * (double precision). The function centers, scales, then calculates the dot
 * product of the pixel with the eigenvector for each output component.
 */
template <typename T>
void 
writePCA(
	std::vector<RasterBandMetaData>& bands,
	std::vector<RasterBandMetaData>& PCABands,
	PCAResult<T>& result,
	GDALDataType type,
	size_t size,
	int height,
	int width)
{
	int bandCount = static_cast<int>(bands.size());
	int nComp = static_cast<int>(PCABands.size());
	std::vector<T> noDataVals(bandCount); 
	std::vector<T *> inputBuffers(bandCount);
	std::vector<T *> PCABandBuffers(nComp);
	std::vector<T *> eigBuffers(nComp);
	T resultNan = std::nan("");
	
	//set no data values and read input bands
	T *p_data = reinterpret_cast<T *>(VSIMalloc3(bandCount, height * width, size);
	for (int b = 0; b < bandCount; b++) {
		noDataVals[b] = static_cast<T>(bands[b].nan);

		inputBuffers[b] = VSIMalloc3(height, width, size);
		bands[b].p_band->RasterIO(
			GF_Read,
			0,
			0,
			width,
			height,
			(void *)((size_t)p_data + b * size),	
			width,
			height,
			type,
			size * bandCount,
			size * bandCount * width
		);	
	}

	//populate PCABandBuffers vector, initialize values to 0, and populate eigBuffers vector
	for (int b = 0; i < nComp; i++) {
		PCABandBuffers[b] = reinterpret_cast<T *>(PCABands[b].p_buffer);
		for (int i = 0; i < height * width; i++) { //this for loop *should* be optimized by compiler
			PCABandBuffers[b][i] = 0;
		}

		eigBuffers[i] = result.eigenvectors[i].data();
	}

	//calculate dot product of eigenvectors and values in a *hopefully* SIMD friendly way
	std::vector<std::vector<T>> multipliers(nComp);
	std::vector<std::vector<T>> adders(nComp);

	for (int c = 0; c < nComp; c++) {
		multipliers[c].resize(bandCount);
		adders[c].resize(bandCount);

		for (int b = 0; b < bandCount; b++) {
			T mean = static_cast<T> result.means[b];
			T stdev = static_cast<T> result.means[b];
			multipliers[c][b] = eigBuffers[c][b] / stdev;
			adders[c][b] = -1 * multipliers[c][b] * mean;
		}
	}
	
	for (int b = 0; b < bandCount; b++) {
		T mean = static_cast<T> result.means[b];
		T stdev = static_cast<T> result.stdevs[b];
		T *buffIn = inputBuffers[b];
		for (int c = 0; c < nComp; c++) {
			T eig = eigBuffers[c][b];
			T mult = eig / stdev;
			T add = -1 * mean * mult;
			buffOut = PCABandBuffers[c];

			for (int i = 0; i < height * width; i++) {
				buffOut[i] += buffIn[i] * mult + add;
				//equivalent to:
				//buffOut[i] += ((buffIn[i] - mean) / stdev) * 	eig
				//
				//which can be decomposed to
				//buffOut[i] += (buffIn[i] - mean) * (eig / stdev)
				//
				//then
				//buffOut[i] += (buffIn[i]) * (eig / stdev) + (-mean * (eig / stdev))
				//
				//setting (eig / stdev) as 'mult' and (-mean * (eig / stdev)) as 'add'
				//limits the number of instructions, and especially division instructions
			}
		}
	}

	//set nan values
	for (int i = 0; i < height * width; i++) {
		for (int b = 0; b < bandCount ; b++) {
			T val = inputBuffers[b][i];
			if (val == noDataVals[b] || std::isnan(val)) {
				for (int c = 0; c < nComp; c++) {
					pcaBandBuffers[c][i] = resultNan;
				}
				break;
			}
		}
	}

	for (int b = 0; b < bandCount; b++) {
		VSIFree(inputBuffers[b]);
	}
}

/**
 * This function is used to write the output principal components to a
 * raster dataset, after the eigenvectors and eigenvalues have already
 * been calculated for the input raster. This function is used in the
 * case where the raster is large, and should be processed in blocks.
 *
 * For each block:
 *
 * First, the input raster band blockss are read into memory using the GDALRasterBand 
 * RasterIO function. Bands are read into memory in a row-wise manor 
 * such that a row indicates a single pixel, and a column indicates a raster band.
 * This means that in between each pixel and the next, a gap must be left for the
 * remaining band values for that pixel index to be written to. This is done
 * using the nPixelSpace, and nLineSpace arguments of RasterIO.
 *
 * Second, the pixels are iterated over. If any value in any band is a no-data
 * value, then nan is written. If none of the values for any band are nodata,
 * then the output pca value is calculated. A different function is called
 * depending on whether the data type is float (single precision), or double 
 * (double precision). The function centers, scales, then calculates the dot
 * product of the pixel with the eigenvector for each output component.
 *
 * Lastly, the output values are written to the output dataset.
 */
template <typename T>
void 
writePCA(
	std::vector<RasterBandMetaData>& bands,
	std::vector<RasterBandMetaData>& PCABands,
	PCAResult<T>& result,
	GDALDataType type,
	size_t size,
	int xBlockSize,
	int yBlockSize,
	int xBlocks,
	int yBlocks)
{
	int bandCount = static_cast<int>(bands.size());
	int nComp = static_cast<int>(PCABands.size());
	std::vector<T *> PCABandBuffers(nComp);
	std::vector<T *> eigBuffers(nComp);
	std::vector<T> noDataVals(bandCount); 
	T resultNan = std::nan("");
	for (int i = 0; i < bandCount; i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}
	for (int i = 0; i < nComp; i++) {
		PCABandBuffers[i] = VSIMalloc3(xBlockSize, yBlockSize, size);
		eigBuffers[i] = reinterpret_cast<void *>(result.eigenvectors[i].data());
	}

	void *p_data = VSIMalloc3(xBlockSize * yBlockSize, size, bandCount);

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
			int xValid, yValid;
			bands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);
		
			//read bands into memory	
			for (int i = 0; i < bandCount; i++) {
				bands[i].p_band->RasterIO(
					GF_Read,
					xBlock * xBlockSize,
					yBlock * yBlockSize,
					xValid,
					yValid,
					(void *)((size_t)p_data + i * size),
					xValid,
					yValid,
					type,
					size * bands.size(),
					size * bands.size() * static_cast<size_t>(xBlockSize)
				);				
			}

			for (int y = 0; y < yValid; y++) {
				int i = y * xBlockSize;
				for (int x = 0; x < xValid; x++) {
					bool isNan = false;
					for (int  b = 0; b < bandCount; b++) {
						T val = reinterpret_cast<T *>(p_data)[i * bandCount + b];
						isNan = val == noDataVals[b] || std::isnan(val);
						if (isNan) {
							break;
						}
					}

					if (isNan) {
						for (void *& p_buffer : PCABandBuffers) {
							reinterpret_cast<T *>(p_buffer)[i] = resultNan;
						}
					}
					else {
						//depending on type, call single or double precision floating point
						//processing functions which use cblas_sdot and cblas_ddot respectively
						type == GDT_Float32 ?
							processSPPixel(i, bandCount, nComp, p_data, PCABandBuffers, eigBuffers, result.means, result.stdevs, n, incx, incy) :
							processDPPixel(i, bandCount, nComp, p_data, PCABandBuffers, eigBuffers, result.means, result.stdevs, n, incx, incy);
					}

					i++;
				}
			}
			
			//write bands to disk
			for (int b = 0; b < nComp; b++) {
				rasterBandIO(
					PCABands[b],
					PCABandBuffers[b],
					xBlockSize,
					yBlockSize,
					xBlock,
					yBlock,
					xValid,
					yValid,
					false, //read = false
					false //threaded = false
				);
			}		
		}

		if (yBlock % 1000 == 0) {
			std::cout << "completed y block " << yBlock << "/" << yBlocks << std::endl;
		}
	}

	VSIFree(p_data);
	for (int i = 0; i < nComp; i++) {
		VSIFree(PCABandBuffers[i]);
	}
}

/**
 * This function conducts principal component analysis on the input raster,
 * writing output bands to a new GDALRasterWrapper, and returning the
 * eigenvectors and eigenvalues calculated for each raster band. The output
 * values are both centered and scaled before being projected onto the pca
 * eigenvectors.
 *
 * First, depending on whether the raster is large (should be processed in
 * blocks) or not, and whether an output filename is given, an output
 * dataset is created to store the output results. In the case of a small
 * raster without a given filename, an in-memory raster is created. In the
 * case of a large raster without a given filename, a VRT dataset is created
 * where each VRT band is a GTiff raster. When a filename is created, the 
 * driver which corresponds to that filename is used.
 *
 * Then, the calculatePCA() function is called, with specific template
 * parameters depending on the data type, and a specific function overload
 * depending on whether the raster should be processed by blocks. This
 * function calculates the principal component eigenvectors, eigenvalues,
 * mean per band, and standard deviation per band. The writePCA() function
 * is then called (again with specific template and overload) to center,
 * scale, and project the input raster values to output pca bands which
 * are written to the output dataset.
 *
 * Finally, a GDALRasterWrapper is created using the output dataset,
 * and returned in a tuple alongside the eigenvectors and eigenvalues. 
 */
std::tuple<GDALRasterWrapper *, std::vector<std::vector<double>>, std::vector<double>>
pca(
	GDALRasterWrapper *p_raster,
	int nComp,
	bool largeRaster,
	std::string tempFolder,
	std::string filename,
	std::map<std::string, std::string> driverOptions)
{
	GDALAllRegister();

	int bandCount = p_raster->getBandCount();
	int height = p_raster->getHeight();
	int width = p_raster->getWidth();
	double *geotransform = p_raster->getGeotransform();
	std::string projection = std::string(p_raster->getDataset()->GetProjectionRef());
	
	bool isMEMDataset = !largeRaster && filename == "";
	bool isVRTDataset = largeRaster && filename == "";
	GDALDataset *p_dataset = nullptr;
	
	std::vector<RasterBandMetaData> bands(bandCount);
	std::vector<RasterBandMetaData> pcaBands(nComp);
	std::vector<VRTBandDatasetInfo> VRTBandInfo;

	int xBlockSize, yBlockSize;
	p_raster->getRasterBand(0)->GetBlockSize(&xBlockSize, &yBlockSize);

	GDALDataType type = GDT_Float32;
	size_t size = sizeof(float);
	for (int i = 0; i < p_raster->getBandCount(); i++) {
		bands[i].p_band = p_raster->getRasterBand(i);
		bands[i].nan = bands[i].p_band->GetNoDataValue();

		if (p_raster->getRasterBandType(i) == GDT_Float64) {
			type = GDT_Float64;
			size = sizeof(double);
		}
	}

	if (isMEMDataset) {
		p_dataset = createVirtualDataset("MEM", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			pcaBands[i].nan = std::nan("");
			addBandToMEMDataset(p_dataset, pcaBands[i]);
		}
	}
	else if (isVRTDataset){
		p_dataset = createVirtualDataset("VRT", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			pcaBands[i].nan = std::nan("");
			createVRTBandDataset(p_dataset, pcaBands[i], tempFolder, pcaBands[i].name, VRTBandInfo, driverOptions);
		}
	}
	else {
		std::filesystem::path filepath = filename;
		std::string extension = filepath.extension().string();
		std::string driver;

		if (extension == ".tif") {
			driver = "GTiff";
		}
		else {
			throw std::runtime_error("sgs only supports .tif files right now.");
		}
		
		bool useTiles = xBlockSize != width &&
				yBlockSize != height;

		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
			pcaBands[i].nan = std::nan("");
			pcaBands[i].p_buffer = !largeRaster ? VSIMalloc3(height, width, size) : nullptr;
		
			if (useTiles) {
				pcaBands[i].xBlockSize = xBlockSize;
				pcaBands[i].yBlockSize = yBlockSize;
			}
		}

		p_dataset = createDataset(
			filename,
			driver,
			width,
			height, 
			geotransform,
			projection,
			pcaBands.data(),
			pcaBands.size(),
			useTiles,
			driverOptions
		);
	}

	//get block size
	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	//these vectors are not used for calculation, but set and returned to the Python side of the application for reference
	std::vector<std::vector<double>> eigenvectors; 
	std::vector<double> eigenvalues;

	//calculate PCA eigenvectors (and eigenvalues), and write values to PCA bands
	switch(type) {
		case GDT_Float32: {
			PCAResult<float> result;
			if (largeRaster) {
				result = calculatePCA<float>(bands, type, size, xBlockSize, yBlockSize, xBlocks, yBlocks, nComp);
				writePCA<float>(bands, pcaBands, result, type, size, xBlockSize, yBlockSize, xBlocks, yBlocks);
			}
			else {
				result = calculatePCA<float>(bands, type, size, width, height, nComp);
				writePCA<float>(bands, pcaBands, result, type, size, height, width);
			}

			eigenvectors.resize(result.eigenvectors.size());
			for (size_t i = 0; i < result.eigenvectors.size(); i++) {
				eigenvectors[i].resize(result.eigenvectors[i].size());
				for (size_t j = 0; j < result.eigenvectors[i].size(); j++) {
					eigenvectors[i][j] = static_cast<double>(result.eigenvectors[i][j]);
				}
			}

			eigenvalues.resize(result.eigenvalues.size());
			for (size_t i = 0; i < result.eigenvalues.size(); i++) {
				eigenvalues[i] = static_cast<double>(result.eigenvalues[i]);
			}
			
			break;
		}
		case GDT_Float64: {
			PCAResult<double> result;
			if (largeRaster) {
				result = calculatePCA<double>(bands, type, size, xBlockSize, yBlockSize, xBlocks, yBlocks, nComp);
				writePCA<double>(bands, pcaBands, result, type, size, xBlockSize, yBlockSize, xBlocks, yBlocks);
			}
			else {
				result = calculatePCA<double>(bands, type, size, width, height, nComp);
				writePCA<double>(bands, pcaBands, result, type, size, height, width);
			}

			eigenvectors = result.eigenvectors;
			eigenvalues = result.eigenvalues;
			
			break;
		}
		default:
			throw std::runtime_error("should not be here! GDALDataType should be one of Float32/Float64!");
	}
	
	//write to file if data is still in memory (the case where largeRaster is false)
	if (!largeRaster && filename != "") {
		CPLErr err;
		for (int b = 0; b < nComp; b++) {
			err = pcaBands[b].p_band->RasterIO(
				GF_Write,
				0,
				0,
				width,
				height, 
				pcaBands[b].p_buffer,
				width,
				height,
				pcaBands[b].type,
				0,
				0
			);
			if (err) {
				throw std::runtime_error("error writing band to file.");
			}
		}
	}

	if (isVRTDataset) {
		for (int b = 0; b < bandCount; b++) {
			GDALClose(VRTBandInfo[b].p_dataset);
			addBandToVRTDataset(p_dataset, pcaBands[b], VRTBandInfo[b]);
		}
	}

	std::vector<void *> buffers(bandCount);
	if (!largeRaster) {
		for (int b = 0; b < bandCount; b++) {
			buffers[b] = pcaBands[b].p_buffer;
		}
	}

	GDALRasterWrapper *p_outrast = largeRaster ?
		new GDALRasterWrapper(p_dataset) :
		new GDALRasterWrapper(p_dataset, buffers);

	return {p_outrast, eigenvectors, eigenvalues};
}

PYBIND11_MODULE(pca, m) {
	m.def("pca_cpp", &pca);
}
