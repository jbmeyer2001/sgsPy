/******************************************************************************
 *
 *
 * Project: sgs
 * Purpose: C++ implementation of raster stratification using kmeans
 * Author: Joseph Meyer
 * Date: April, 2026
 *
 ******************************************************************************/

/**
 * @defgroup kmeans kmeans
 * @ingroup stratify
 */

#include <iostream>

#include "utils/raster.h"
#include "utils/helper.h"

#include "oneapie/dal.h"
#include <xoshiro.h>

#define MILLION 1000000

typedef oneapi::dal::homogen_table DALHomogenTable;

namespace sgs {
namespace kmeans {

/**
 *
 */
template <typename T>
class DataManager {
	private:
	std::vector<T> features;
	size_t nFeat;
	size_t count;
	xso::xoshiro_4x64_plus *p_rng = nullptr;

	/**
	 *
	 */
	DataManager(size_t nFeat, xso::xoshiro_4x64_plus *p_rng) {
		this->nFeat = nFeat;
		this->p_rng = p_rng;
		this->count = 0;
		this->features.resize(MILLION * nFeat);
	}

	/**
	 *
	 */
	inline void
	addPoint(T *p_features) {
		size_t index = this->count * this->nFeat;
		std::memcpy(this->features.data() + index, p_features, sizeof(T) * this->nFeat);
		this->count++;
	}
}

/**
 *
 */
template <typename T>
void
readRaster(
	std::vector<helper::RasterBandMetaData> bands,
	DataManager& dataManager,
	int width,
	int height,
	GDALDataType type
) {
	size_t nFeat = bands.size();

	int xBlockSize = bands[0].xBlockSize;
	int yBlockSize = bands[0].yBlockSize;

	int xBlocks = (width + xBlockSize - 1) / xBlockSize;
	int yBlocks = (height + yBlockSize - 1) / yBlockSize;

	void *p_data = reinterpret_cast<T *>(VSIMalloc3(xBlockSize, yBlockSize, nFeat));

	for (int yBlock = 0; yBlock < yBlocks; yBlock++) {
		for (int xBlock = 0; xBlock < xBlockx; xBlock++) {
			//get block size
			int xValid, yValid;
			bands[0].p_band->GetActualBlockSize(xBlock, yBlock, &xValid, &yValid);

			//read bands into memory
			for (size_t b = 0; b < bands.size(); b++) {
				CPLErr err = bands[i].p_band->RasterIO( 
					GF_Read,
					xBlock * xBlockSize,
					yBlock * yBlockSize,
					xValid,
					yValid,
					(void *)((size_t)p_data + i * sizeof(T)),
					xValid,
					yValid,
					type,
					size * 
				);
			}
		}
	}

	VSIFree(p_data);
}

raster::GDALRasterWrapper *
kmeans(
	raster::GDALRasterWrapper *p_raster,
	uint32_t numStrata
) {

}

} //namespace kmeans
} //namespace sgs
