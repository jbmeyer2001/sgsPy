# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratified random sampling (srs)
#  Author: Joseph Meyer
#  Date: September, 2025
#
# ******************************************************************************

##
# @defgroup user_strat strat
# @ingroup user_sample

import os
import sys
import site
import tempfile
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgspy.utils import(
    SpatialRaster,
    SpatialVector,
    StratRasterBandMetadata,
    plot,
)

#ensure _sgs binary can be found
site_packages = list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0]
sys.path.append(os.path.join(site_packages, "sgspy"))
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from _sgs import strat_cpp

##
# @ingroup user_strat
# This function conducts stratified sampling using the stratified
# raster given. There are two methods employed to determine which
# pixels to sample:
#  - The 'random' method randomly selects pixels within a given strata.
#  - The 'Queinnec' method prioritizes pixels which are surrounded by other pixels of the same strata.
# The 'wrow' and 'wcol' parameters determine the size of the surrounding area required for a pixel to be
# prioritized, and must be one of the following integers: 3, 5, 7.
# 
# The desired number of samples is given by num_samples.
#
# **IMPORTANT**
# the num_strata argument is required only if the strat raster given is not the return value of a sgspy
# stratification function. If the strat raster given is the return value of an sgspy stratification 
# function, then this value is automatically stored and used. The num_strata argument, if required, should
# be set to the value of the largest strata + 1. For example if the strata are [0, 1, 2, 3, 4] then num_strata should
# be 5. If the strata are [1, 2, 4] then num_strata should still be 5. If the strata are [0, 1, 2, 3] then
# num_strata should be 4.
#
# The allocation parameter specifies the proportion of total samples will be distributed between
# each strata. The 'prop' method is the default, and attempts to allocate the samples proportionally according to 
# their prevalence in the overall raster. The 'equal' method attempts to distribute the samples equally among strata.
# the 'manual' method requires that the weights parameter be given, and attempts to allocate according to the
# proportions given in the weights parameter. In the case where 'optim' allocation is used, an additional raster must be passed
# to the mrast parameter, and if that raster contains more than 1 band the mrast_band
# parameter must be given specifying which band. The optim method is specified by Gregoire and Valentine,
# and optimizes the desired proportions based on the proportion of each strata AND the within-strata
# variance in the specified raster band. https://doi.org/10.1201/9780203498880 Section 5.4.4.
# 
# The 'existing' parameter, if passed, must be a SpatialVector of type Point or MultiPoint. 
# These points specify samples within an already-existing network. The SpatialVector may
# only have one layer. If the force parameter is set to True, every pixel in the existing
# sample will be added no matter what. if the force parameter is false, then the existing
# samples will be prioritized over other pixels in the same strata.
# 
# The 'access' parameter, if passed, must be a SpatialVector of type LineString or MultiLineString.
# buff_outer specifies the buffer distance around the geometry which
# is allowed to be included in the sampling, buff_inner specifies the
# geometry which is not allowed to be included in the sampling. buff_outer
# must be larger than buff_inner. For a multi-layer vector, layer_name
# must be specified.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster("raster.tif") @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5) @n
# samples = sgspy.sample.strat(srast, num_samples=200) #uses Queinnec method with proportional allocation by default
# 
# srast = sgs.SpatialRaster("srast.tif") #srast not result of sgspy stratification function, num_strata arg necessary in sgspy.sample.strat
# samples = sgspy.sample.strat(srast, num_strata=5, num_samples=200)
#
# rast = sgspy.SpatialRaster("raster.tif") @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5) @n
# samples = sgspy.sample.strat(srast, num_samples=200, method="random", mindist=200, plot=True, filename="samples.shp")
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5) @n
# samples = sgspy.sample.strat(srast, num_samples=200, method="Queinnec", allocation="optim", mrast=rast)
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5) @n
# samples = sgspy.sample.strat(rast, num_samples=200, allocation="manual", weights=[0.1, 0.1, 0.2, 0.2, 0.4])
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# access = sgspy.SpatialVector("access_network.shp") @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5) @n
# samples = sgspy.sample.strat(rast, num_samples=200, allocation="equal", access=access, buff_inner=100, buff_outer=300)
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# existng = sgspy.SpatialVector("existing_samples.shp") @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5) @n
# samples = sgspy.sample.strat(rast, num_samples=200, allocation="prop", existing=existing, force=True)
# 
# Parameters
# --------------------
# strat_rast : SpatialRaster
#     the raster to sample
# band : Optional[int | str] @n
#     the band within the strat_rast to use, required if strat_rast has more than 1 band, either a 0-indexed int value or the name of the band @n @n
# num_samples : int @n
#     the desired number of samples @n @n
# num_strata : Optional[int] @n
#     the value of the largest stratification in the strat_rast + 1 @n @n
# wrow : int @n
#     the number of rows to be considered in the focal window for the 'Queinnec' method @n @n
# wcol : int @n
#     the number of columns to be considered in the focal window for the 'Queinnec' method @n @n
# allocation : str @n
#     the allocation method to determine the number of samples per strata. One of 'prop', 'equal', 'optim', or 'manual' @n @n
# weights : list[float] @n
#     the allocation percentages of each strata if the allocation method is 'manual' @n @n
# mrast : SpatialRaster @n
#     the raster used to calculate 'optim' allocation if 'optim' allocation is used @n @n
# mrast_band : str | int @n
#     specifies the band within mrast to use @n @n
# method : str @n
#     the sampling method, either 'random', or 'Queinnec' @n @n
# mindist : float @n
#     the minimum distance allowed between sample points @n @n
# existing : SpatialVector @n
#     a vector of Point or Multipoint which are part of a pre-existing sample network @n @n
# force : bool @n
#     whether to automatically include all points in the existing network or not @n @n
# access : SpatialVector @n
#     a vector of LineString or MultiLineString geometries to sample near to @n @n
# layer_name : str @n @n
#     the layer within 'access' to use for access buffering @n @n
# buff_inner : float @n
#     the inner buffer around the access LineStrings, where samples should not occur @n @n
# buff_outer : float @n
#     the outer buffer around the access LineStrings, The area in which samples must occur @n @n
# plot : bool @n
#     whether or not to plot the output samples @n @n
# filename : str @n
#     the output filename to write to if desired @n @n
# 
# 
# Returns
# --------------------
# a SpatialVector object containing point geometries of sample locations
def strat(
    strat_rast: SpatialRaster,
    num_samples: int,
    num_strata: Optional[int] = None,
    wrow: int = 3,
    wcol: int = 3,
    band: Optional[int | str] = None,
    allocation: str = "prop",
    weights: Optional[list[float]] = None,
    mrast: Optional[SpatialRaster] = None,
    mrast_band: Optional[int | str] = None,
    method: str = "Queinnec",
    mindist: Optional[int | float] = None,
    existing: Optional[SpatialVector] = None,
    force: bool = False,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    plot: bool = False,
    filename: str = "",
    ):

    if type(strat_rast) is not SpatialRaster:
        raise TypeError("'strat_rast' parameter must be of type sgspy.SpatialRaster.")

    if band is not None and type(band) not in [int, str]:
        raise TypeError("'band' parameter, if given, must be of type int or str.")

    if type(num_samples) is not int:
        raise TypeError("'num_samples' parameter must be of type int.")

    if num_strata is not None and type(num_strata) is not int:
        raise TypeError("'num_strata' parameter, if given, must be of type int.")

    if not strat_rast.is_strat_rast and num_strata is None:
        raise ValueError("'if 'strat_rast' parameter is not the return value of an sgspy.stratify function, 'num_strata' is a required parameter.")

    if type(wrow) is not int:
        raise TypeError("'wrow' parameter must be of type int.")

    if type(wcol) is not int:
        raise TypeError("'wcol' parameter must be of type int.")

    if type(allocation) is not str:
        raise TypeError("'allocation' parameter must be of type str.")

    if weights is not None and type(weights) is not list:
        raise TypeError("'weights' parameter, if given, must be a list of float values.")

    if mrast is not None and type(mrast) is not SpatialRaster:
        raise TypeError("'mrast' parameter, if given, must be of type sgspy.SpatialRaster.")

    if mrast_band is not None and type(mrast_band) not in [int, str]:
        raise TypeError("'mrast_band' parameter, if given, must be of type int or str.")

    if type(method) is not str:
        raise TypeError("'method' parameter must be of type str.")

    if mindist is not None and type(mindist) not in [int, float]:
        raise TypeError("'mindist' parameter must be of type int or float.")

    if existing is not None and type(existing) is not SpatialVector:
        raise TypeError("'existing' parameter must be of type sgspy.SpatialVector.")

    if type(force) is not bool:
        raise TypeError("'force' parameter must be of type bool.")

    if access is not None and type(access) is not SpatialVector:
        raise TypeError("'access' parameter must be of type sgspy.SpatialVector.")

    if layer_name is not None and type(layer_name) is not str:
        raise TypeError("'layer_name' parameter must be of type str.")

    if buff_inner is not None and type(buff_inner) not in [int, float]:
        raise TypeError("'buff_inner' parameter must be of type int or float.")

    if buff_outer is not None and type(buff_outer) not in [int, float]:
        raise TypeError("'buff_outer' parameter must be of type int or float.")

    if type(plot) is not bool:
        raise TypeError("'plot' paramter must be of type bool.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if strat_rast.closed:
        raise RuntimeError("the C++ object which the strat_rast object wraps has been cleaned up and closed.")

    if mrast is not None and mrast.closed:
        raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    if band is None:
        if len(strat_rast.bands) > 1:
            raise ValueError("'band' parameter must be given if there is more than 1 band in the strat_raster")
        band = 0
    elif type(band) is str:
        if band not in strat_rast.bands:
            msg = "band " + str(band) + " not in given raster."
            raise ValueError(msg);

        band = strat_rast.get_band_index(band) #get band as 0-indexed integer
    else: #type(band) is int
        if band >= len(strat_rast.bands):
            msg = "0-indexed band of " + str(band) + " given, but raster only has " + str(len(raster.bands)) + " bands."
            raise ValueError(msg)

    if strat_rast.is_strat_rast:
        num_strata = strat_rast.srast_metadata_info[band].get_num_strata() 

    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

    if method not in ["random", "Queinnec"]:
        raise ValueError("method must be either 'random' or 'Queinnec'.")

    if allocation not in ["prop", "optim", "equal", "manual"]:
        raise ValueError("method must be one of 'prop', 'optim', 'equal', or 'manual'.")

    if allocation == "manual":
        if weights is None:
            raise ValueError("for manual allocation, weights must be given.")

        if np.sum(weights) != 1:
            raise ValueError("weights must sum to 1.")

        if len(weights) != num_strata:
            raise ValueError("length of 'weights' must be the same as the number of strata, which is {}".format(num_strata))

    if allocation == "optim":
        if not mrast:
            raise ValueError("the 'mrast' parameter must be provided if a SpatialRaster if allocation is 'optim'.")

        if mrast_band is None:
            if len(mrast.bands) != 1:
                raise ValueError("the 'mrast_band' parameter must be given if the 'mrast' SpatialRaster contains more than 1 band.")

            mrast_band = 1
        elif type(mrast_band) is str:
            if band not in mrast.bands:
                msg = "band " + str(band) + " not in mraster."
                raise ValueError(msg)

            mrast_band = mrast.get_band_index(mrast_band)
        else: #type(band) is int
            if (mrast_band >= len(mrast.bands)):
                msg = "0-indexed band of " + str(band) + "given, but raster only has " + str(len(mrast.bands)) + " bands."
                raise ValueError(msg)

        mrast_cpp_raster = mrast.cpp_raster
    else:
        mrast_cpp_raster = None
        mrast_band = -1

    if wrow not in [3, 5, 7]:
        raise ValueError("wrow must be one of 3, 5, 7.")

    if wcol not in [3, 5, 7]:
        raise ValueError("wcol must be one of 3, 5, 7.")

    if (access):
        if layer_name is None:
            if len(access.layers) > 1:
                raise ValueError("if there are multiple layers in the access vector, layer_name must be defined.")
            layer_name = access.layers[0]

        if layer_name not in access.layers:
            raise ValueError("layer specified by 'layer_name' does not exist in the access vector")

        if buff_inner is None or buff_inner < 0:
            buff_inner = 0

        if buff_outer is None or buff_outer <= 0:
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

    if allocation != "manual":
        weights = []

    if mindist is None:
        mindist = 0

    if mindist < 0:
        raise ValueError("mindist must be greater than or equal to 0")

    temp_dir = strat_rast.cpp_raster.get_temp_dir()
    if temp_dir == "":
        temp_dir = tempfile.mkdtemp()
        strat_rast.cpp_raster.set_temp_dir(temp_dir)

    [sample_coordinates, samples, num_points] = strat_cpp(
        strat_rast.cpp_raster,
        band,
        num_samples,
        num_strata,
        allocation,
        weights,
        mrast_cpp_raster,
        mrast_band,
        method,
        wrow,
        wcol,
        mindist,
        existing_vector,
        force,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        plot,
        filename,
        temp_dir
    )

    if num_points < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, num_points))
    
    #plot new vector if requested
    if plot:
        try:
            fig, ax = plt.subplots()
            strat_rast.plot(ax, band=strat_rast.bands[band])
            title = "samples on " + strat_rast.bands[band]

            if access:
                access.plot('LineString', ax)
                title += " with access"

            ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
            ax.set_title(label=title)
            plt.show()
        except Exception as e:
            print("unable to plot output: " + str(e))

    return SpatialVector(samples)
