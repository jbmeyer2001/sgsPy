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
#include <mkl.h>

typedef oneapi::dal::homogen_table	DALHomogenTable;

/**
 *
 */
template <typename T>
struct PCAResult {
	std::vector<std::vector<T>> eigenvectors;
	std::vector<T> eigenvalues;
	std::vector<double> means;
	std::vector<double> stdevs;
};

/**
 *
 */
//https://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
struct Variance {
	int64_t k;
	double M = 0;
	double S = 0;
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
		double variance = S / static_cast<double>(k - 1);
		return std::sqrt(variance);
	}
};

/**
 *
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
 *
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

	for (int xBlock = 0; xBlock < xBlocks; xBlock++) {
		for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
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
			int nFeatures;
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
 *
 */
inline void 
processSPPixel(
	int i,
	int bandCount,
	int nComp,
	void *p_data,
	std::vector<void *>& PCABandBuffers,
	std::vector<void *>& eigBuffers,
	std::vector<double>& means,
	std::vector<double>& stdevs,
	MKL_INT n,
	MKL_INT incx,
	MKL_INT incy)
{
	//get features array
	float *p_features = reinterpret_cast<float *>(p_data) + i * bandCount;
	
	//scale and center features
	for (int b = 0; b < bandCount; b++) {
		p_features[b] = (p_features[b] - static_cast<float>(means[b])) / static_cast<float>(stdevs[b]);
	}

	//use blas to calculate projection (dot product) for output raster
	for (int b = 0; b < nComp; b++) {
		reinterpret_cast<float *>(PCABandBuffers[b])[i] = cblas_sdot(
			n,
			p_features,
			incx,
			reinterpret_cast<float *>(eigBuffers[b]),
			incy
		);
	}
}

/**
 *
 */
inline void 
processDPPixel(
	int i,
	int bandCount,
	int nComp,
	void *p_data,
	std::vector<void *>& PCABandBuffers,
	std::vector<void *>& eigBuffers,
	std::vector<double>& means,
	std::vector<double>& stdevs,
	MKL_INT n,
	MKL_INT incx,
	MKL_INT incy)
{
	//get features array
	double *p_features = reinterpret_cast<double *>(p_data) + i * bandCount;
	
	//scale and center features
	for (int b = 0; b < bandCount; b++) {
		p_features[b] = (p_features[b] - means[b]) / stdevs[b];
	}

	//use blas to calculate projection (dot product) for output raster
	for (int b = 0; b < nComp; b++) {
		reinterpret_cast<double *>(PCABandBuffers[b])[i] = cblas_ddot(
			n,
			p_features,
			incx,
			reinterpret_cast<double *>(eigBuffers[b]),
			incy
		);
	}
}


/**
 *
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
	void *p_data = VSIMalloc3(width * height, bandCount, size);
	std::vector<void *> PCABandBuffers(nComp);
	std::vector<void *> eigBuffers(nComp);
	std::vector<T> noDataVals(bandCount); 
	T resultNan = -1;
	for (int i = 0; i < bandCount; i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}
	for (int i = 0; i < nComp; i++) {
		PCABandBuffers[i] = reinterpret_cast<void *>(PCABands[i].p_buffer);
		eigBuffers[i] = reinterpret_cast<void *>(result.eigenvectors[i].data());
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

	MKL_INT n = result.eigenvectors.size();
	MKL_INT incx = 1;
	MKL_INT incy = 1;

	//process chunk of data
	for (int i = 0; i < height * width; i++) {
		bool isNan;
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
	}
}

/**
 *
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
	std::vector<void *> PCABandBuffers(nComp);
	std::vector<void *> eigBuffers(nComp);
	std::vector<T> noDataVals(bandCount); 
	T resultNan = -1;
	for (int i = 0; i < bandCount; i++) {
		noDataVals[i] = static_cast<T>(bands[i].nan);
	}
	for (int i = 0; i < nComp; i++) {
		PCABandBuffers[i] = reinterpret_cast<void *>(PCABands[i].p_buffer);
		eigBuffers[i] = reinterpret_cast<void *>(result.eigenvectors[i].data());
	}

	void *p_data = VSIMalloc3(xBlockSize * yBlockSize, size, bandCount);

	MKL_INT n = result.eigenvectors.size();
	MKL_INT incx = 1;
	MKL_INT incy = 1;

	for (int yBlock = 0; yBlock < yBlockSize; yBlock++) {
		for (int xBlock = 0; xBlock < xBlockSize; xBlock++) {
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
					bool isNan;
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
			for (int b = 0; b < bandCount; b++) {
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
	}
}

/**
 *
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
			addBandToMEMDataset(p_dataset, pcaBands[i]);
		}
	}
	else if (isVRTDataset){
		p_dataset = createVirtualDataset("VRT", width, height, geotransform, projection);
	
		for (int i = 0; i < nComp; i++) {
			pcaBands[i].type = type == GDT_Float64 ? GDT_Float64 : GDT_Float32;
			pcaBands[i].size = type == GDT_Float64 ? sizeof(double) : sizeof(float);
			pcaBands[i].name = "comp_" + std::to_string(i + 1);
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
