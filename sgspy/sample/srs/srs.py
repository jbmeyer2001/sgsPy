# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

##
# @defgroup user_srs srs
# @ingroup user_sample

import os
import sys
import tempfile
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgspy.utils import (
    SpatialRaster,
    SpatialVector,
    plot,
)

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from _sgs import srs_cpp

##
# @ingroup user_srs
# This function conducts simple random sampling on the raster given. 
# Sample points are randomly selected from data pixels (can't be nodata).
# All sample points are at least mindist distance away from eachother.
# If unable to get the full number of sample points, a message is printed.
#
# An access vector of LineString or MultiLineString type can be provided.
# buff_outer specifies the buffer distance around the geometry which
# is allowed to be included in the sampling, buff_inner specifies the
# buffer distance around the geometry which is not allowed to be included 
# in the sampling. buff_outer must be larger than buff_inner. For a multi 
# layer vector, layer_name must be specified.
#
# A vector containing existing sample points can be provided. If this is
# the case then all of the points in the existing sample are automatically
# added and random samples are chosen as required until num_samples number
# of samples are chosen.
#
# Examples
# --------------------
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.srs(rast, num_samples=250) 
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.srs(rast, num_samples=250, mindist=100, plot=True, filename="srs_samples.shp") @n
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("access_network.shp") @n
# samples = sgspy.sample.srs(rast, num_samples=200, mindist=100, access=access, buff_outer=300)
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("access_network.shp") @n
# samples = sgspy.sample.srs(rast, num_samples=200, access=access, buff_inner=50, buff_outer=300)
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# existing = sgspy.SpatialVector("existing_samples.shp") @n
# samples = sgspy.sample.srs(rast, num_samples=200, existing=existing)
#
# Parameters
# --------------------
# rast : SpatialRaster @n
#     raster data structure containing the raster to sample @n @n
# num_samples : int @n
#     the target number of samples @n @n
# mindist : float @n
#     the minimum distance each sample point must be from each other @n @n
# existing : SpatialVector @n
#     a vector specifying existing sample points @n @n
# access : SpatialVector @n
#     a vector specifying access network @n @n
# layer_name : str @n
#     the layer within access that is to be used for sampling @n @n
# buff_inner : int | float @n
#     buffer boundary specifying distance from access which CANNOT be sampled @n @n
# buff_outer : int | float @n
#     buffer boundary specifying distance from access which CAN be sampled @n @n
# plot : bool @n
#     whether to plot the samples or not @n @n
# filename : str @n
#     the filename to write to, or '' if file should not be written @n @n
#
#
# Returns
# --------------------
# a SpatialVector object containing point geometries of sample locations
def srs(
    rast: SpatialRaster,
    num_samples: int,
    mindist: [int | float] = 0,
    existing: Optional[SpatialVector] = None,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    plot: bool = False,
    filename: str = ''):
        
    if type(rast) is not SpatialRaster:
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster.")

    if type(num_samples) is not int:
        raise TypeError("'num_samples' parameter must be of type int.")

    if type(mindist) not in [int, float]:
        raise TypeError("'mindist' parameter must be of type int or float.")

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

    if type(plot) is not bool:
        raise TypeError("'plot' parameter must be of type bool.")

    if type(filename) is not str:
        raise TypeError("'filename' paramter must be of type str.")

    if rast.closed:
            raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

    if mindist is None:
        mindist = 0

    if mindist < 0:
        raise ValueError("mindist must be greater than or equal to 0")



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

    temp_dir = rast.cpp_raster.get_temp_dir()
    if temp_dir == "":
        temp_dir = tempfile.mkdtemp()
        rast.cpp_raster.set_temp_dir(temp_dir)

    #call random sampling function
    [sample_coordinates, cpp_vector, num_points] = srs_cpp(
        rast.cpp_raster,
        num_samples,
        mindist,
        existing_vector,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        plot,
        temp_dir,
        filename
    )
    
    if num_points < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, num_points))

    #plot new vector if requested
    if plot:
        try:
            fig, ax = plt.subplots()
            rast.plot(ax, band=rast.bands[0])
            title = "samples on " + rast.bands[0]
            
            if access:
                access.plot('LineString', ax)
                title += " with access"

            ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
            ax.set_title(label=title)
            plt.show()

        except Exception as e:
            print("unable to plot output: " + str(e))

    return SpatialVector(cpp_vector)
