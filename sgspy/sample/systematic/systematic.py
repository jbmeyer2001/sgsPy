# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

##
# @defgroup user_systematic systematic
# @ingroup user_sample

import os
import sys
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgspy.utils import (
    SpatialRaster,
    SpatialVector,
    plot,
)

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from _sgs import systematic_cpp

##
# @ingroup user_systematic
# This function conducts systematic sampling within the extent of
# the raster given. The 'cellsize' parameter specifies the grid size,
# the 'shape' parameter specifies the grid shape, and the 'location' 
# parameter specifies where in the grid a sample should fall into.
# 
# shape can be one of 'square', and 'hexagon'.
# location can be one of 'corners', 'centers', 'random'.
# 
# An access vector of LineString or MultiLineString type can be provided.
# buff_outer specifies the buffer distance around the geometry which is
# allowed to be included in the sampling, buff_inner specifies the geometry
# which is not allowed to be included in the sampling. buff_outer must
# be larger than buff_inner. For a multi-layer vector, layer_name
# must be provided.
# 
# A vector containing existing sample points can be provided. If this is
# the case then all of the points in the existing sample are automatically
# added and random samples are then chosen as required until num_samples 
# number of samples are chosen.
# 
# If the force parameter is True, then the the samples are forced to 
# fall on an index which is NOT a no data value. This may result
# in some grids not being sampled.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.systematic(rast, 500, "hexagon", "centers")
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.systematic(rast, 500, "square", "corners", plot=True, filename="systematic_samples.shp")
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.systematic(rast, 500, "hexagon", "random", force=True)
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("access_network.shp") @n
# samples = sgspy.sample.systematic(rast, 500, "hexagon", "corners", access=access, buff_outer=300)
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("existing_samples.shp") @n
# samples = sgspy.sample.systematic(rast, 500, "hexagon", "corners", existing=existing)
# 
# Parameters
# --------------------
# rast : SpatialRaster @n
#     the raster to be sampled @n @n
# cellsize : float @n
#     the size of the grid cells to be sampled @n @n
# shape : str @n
#     the shape of the grid cells to be sampled @n @n
# location : str @n
#     the location within the grid cell to be sampled @n @n
# existing (optional) : SpatialVector @n
#     a vector specifying existing sample points @n @n
# access (optional) : SpatialVector @n
#     a vector specifying access network @n @n
# layer_name (optional) : str @n
#     the layer within access that is to be used for sampling @n @n
# buff_inner (optional) : int | float @n
#     buffer boundary specifying distance from access which CANNOT be sampled @n @n
# buff_outer (optional) : int | float @n
#     buffer boundary specifying distance from access which CAN be sampled @n @n
# force : bool @n
#     True if samples are not allowed to fall on a nodata pixel @n @n
# plot : bool @n
#     whether or not to plot the resulting samples @n @n
# filename : str @n
#     the filename to write to or "" if not to write @n @n
# 
# Returns
# --------------------
# a SpatialVector object containing point geometries of sample locations
def systematic(
    rast: SpatialRaster,
    cellsize: int | float,
    shape: str = "square",
    location: str = "centers",
    existing: Optional[SpatialVector] = None,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    force: bool = False,
    plot: bool = False,
    filename: str = ""):
        
    if type(rast) is not SpatialRaster:
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster.")

    if type(cellsize) not in [int, float]:
        raise TypeError("'cellsize' parameter must be of type int or float.")

    if type(shape) is not str:
        raise TypeError("'shape' paramter must be of type str.")

    if type(location) is not str:
        raise TypeError("'location' parameter must be of type str.")

    if existing is not None and type(existing) is not SpatialVector:
        raise TypeError("'existing' parameter, if given, must be of type sgspy.SpatialVector.")

    if access is not None and type(access) is not SpatialVector:
        raise TypeError("'access' parameter, if given, must be of type sgspy.SpatialVector.")

    if layer_name is not None and type(layer_name) is not str:
        raise TypeError("'layer_name' parameter, if given, must be of type str.")

    if buff_inner is not None and type(buff_inner) not in [int, float]:
        raise TypeError("'buff_inner' parameter, if given, must be of type int or float.")

    if buff_outer is not None and type(buff_outer) not in [int, float]:
        raise TypeError("'buff_outer' parameter, if given, must be of type int or float.")

    if type(force) is not bool:
        raise TypeError("'force' parameter must be of type bool.")

    if type(plot) is not bool:
        raise TypeError("'plot' parameter must be of type bool.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if rast.closed:
        raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    if cellsize <= 0:
        raise ValueError("cellsize must be greater than 0")

    if shape not in ["square", "hexagon"]:
        raise ValueError("shape parameter must be one of 'square', 'hexagon'")

    if location not in ["centers", "corners", "random"]:
        raise ValueError("location parameter must be one of 'centers', 'corners', 'random'")

    if (access):
        if layer_name is None:
            if len(access.layers) > 1:
                raise ValueError("if there are multiple layers in the access vector, layer_name parameter must be passed.")
            layer_name = access.layers[0]

        if layer_name not in access.layers:
            raise ValueError("layer specified by 'layer_name' does not exist in the access vector")

        if buff_inner is None or buff_inner < 0:
            buff_inner = 0

        if buff_outer is None or buff_outer < 0:
            raise ValueError("if an access vector is given, buff_outer must be a float greater than 0.")

        if buff_inner >= buff_outer:
            raise ValueError("buff_outer must be greater than buff_inner")

        access_vector = access.cpp_vector
    else:
        access_vector = None
        layer_name = ""
        buff_inner = -1
        buff_outer = -1

    if (existing):
        existing_vector = existing.cpp_vector
    else:
        existing_vector = None

    [samples, points, grid] = systematic_cpp(
        rast.cpp_raster,
        cellsize,
        shape,
        location,
        existing_vector,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        force,
        plot,
        filename
    )

    #plot new vector if requested
    if plot:
        fig, ax = plt.subplots()
        ax.set_xlim([rast.xmin, rast.xmax])
        ax.set_ylim([rast.ymin, rast.ymax])
        rast.plot(ax, band=rast.bands[0])
        title="samples on " + rast.bands[0]
        
        #plot grid
        for shape in grid:
            ax.plot(shape[0], shape[1], '-k')

        #plot sample points
        ax.plot(points[0], points[1], '.r')
        ax.set_title(label=title)
        plt.show()

    return SpatialVector(samples)
