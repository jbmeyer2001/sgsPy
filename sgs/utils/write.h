/******************************************************************************
 *
 * Project: sgs
 * Purpose: Write sample points to vector file
 * Author: Joseph Meyer
 * Date: June, 2025
 *
 ******************************************************************************/

#pragma once

#include <gdal_priv.h>

/**
 * Write function for writing sample points to a vector file.
 *
 * @param std::vector<OGRPoint> the points to write.
 * @param std::string the filename to write to
 */
void write(std::vector<OGRPoint> points, std::string filename);
