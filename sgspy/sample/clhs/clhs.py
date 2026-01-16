# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

##
# @defgroup user_clhs clhs
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
from _sgs import clhs_cpp

## 
# @ingroup user_clhs
# This function conducts Conditioned Latin Hypercube Sampling, see the following article for an
# in depth description of the method itself:
# 
# Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method
# for sampling in the presence of ancillary information. Computers and Geosciences, 32:1378-1388.
#
# The number of output samples is decided by the 'num_samples' parameter. The 'iterations' parameter
# indicates the number of iterations the simulated annealing portion of the clhs algorithm will undertake
# in the case where a perfect latin hypercube is not found. A higher number of iterations may result in
# a more representative sample, although the standard value recommended by Misany and McBratney is 10000.
#
# The access parameter may be given to restrict the areas where sampling may occur. The algorithm will still
# attempt to find a latin hypercube representative across the entire feature space, not just the accessible
# pixels. The access vector may contain geometries of type LineString or MultiLineString. buff_outer specifies
# the buffer distance around the geometry which is allowed to be included in the sampling. buff_inner specifies
# the buffer distance around the geometry which is not allwoed to be included in the sampling. buff_outer must
# be larger than buff_inner. For a multi layer vector, layer_name must be specified.
#
# The output is an object of type sgspy.SpatialVector which contains the chosen sample points.
#
# Examples
# --------------------
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.clhs(rast, num_samples=250)
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# samples = sgspy.sample.clhs(rast, num_samples=250, plot=True, filename="clhs_samples.shp")
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("access_network.shp") @n
# samples = sgspy.sample.clhs(rast, num_samples=200, access=access, buff_outer=300)
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("access_network.shp") @n
# samples = sgspy.sample.clhs(rast, num_samples=200, access=access, buff_inner=50, buff_outer=300)
#
# Parameters
# --------------------
# rast : SpatialRaster @n
#     raster data structure containing input raster bands @n @n
# num_samples : int @n
#     the target number of samples @n @n
# iterations : int @n
#     the number of iterations in the clhs algorithms @n @n
# access : SpatialVector @n
#     a vector specifying an access network @n @n
# layer_name : str @n
#     the layer within the access network which will be used for sampling @n @n
# buff_inner : int | float @n
#     buffer boundary specifying distance from access geometries which CANNOT be sampled @n @n
# buff_outer : int | float @n
#     buffer boundary specifying distance from access geometries which CAN be sampled @n @n
# plot : bool @n
#     whether to plot the output samples or not @n @n
# filename : str @n
#     the filename to write to, or '' if file should not be written @n @n
#
# Returns
# --------------------
# a SpatialVector object containing point geometries of sample locations
def clhs(
    rast: SpatialRaster,
    num_samples: int,
    iterations: int = 10000,
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

    if type(iterations) is not int:
        raise TypeError("'iterations' parameter must be of type int.")

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
        raise TypeError("'filename' parameter must be of type str.")

    if rast.closed:
            raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

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

    temp_dir = rast.cpp_raster.get_temp_dir()
    if temp_dir == "":
        temp_dir = tempfile.mkdtemp()
        rast.cpp_raster.set_temp_dir(temp_dir)

    [sample_coordinates, cpp_vector] = clhs_cpp(
        rast.cpp_raster,
        num_samples,
        iterations,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        plot,
        temp_dir,
        filename
    )

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
