# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratified random sampling (srs)
#  Author: Joseph Meyer
#  Date: September, 2025
#
# ******************************************************************************

import tempfile
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgs.utils import(
    SpatialRaster,
    SpatialVector,
    plot,
)

from strat import strat_cpp

def strat(
    strat_rast: SpatialRaster, #TODO add band name for strat rast
    band: int | str,
    num_samples: int,
    num_strata: int,
    wrow: int = 3,
    wcol: int = 3,
    allocation: str = "prop",
    weights: Optional[list[float]] = None,
    mrast: Optional[SpatialRaster] = None,
    mrast_band: Optional[int | str] = None,
    method: str = "Queinnec",
    mindist: Optional[float] = None,
    existing: Optional[SpatialVector] = None,
    force: bool = False,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[float] = None,
    buff_outer: Optional[float] = None,
    plot: bool = False,
    filename: str = "",
    ):
    """
    This function conducts stratified sampling using the stratified
    raster given. There are two methods employed to determine which
    pixels to sample:
     - The 'random' method randomly selects pixels
    within a given strata.
     - The 'Queinnec' method first selects pixels which are surrounded
    by pixels of the same strata, the focal window, which is defined by 
    the wrow and wcol parameters.

    The number of total samples is given by num_samples. The allocation
    of samples per strata is calculated given the distribution of pixels
    in each strata, and the allocation method specified by the allocation parameter.

    Parameters
    --------------------
    strat_rast : SpatialRaster
        the raster to sample
    band : int | str
        the band within the strat_rast to use, either a 0-indexed int value or the name of the band
    num_samples : int
        the desired number of samples
    num_strata : Int
        the number of strata in the strat_rast. If this number is incorrect it may 
        undefined behavior in the C++ code which determines sample locations.
    wrow : int
        the number of rows to be considered in the focal window for the 'Queinnec' method
    wcol : int
        the number of columns to be considered in the focal window for the 'Queinnec' method
    allocation : str
        the allocation method to determine the number of samples per strata. One of 'prop', 'equal', 'optim', or 'manual'
    method: str
        the sampling method, either 'random', or 'Queinnec'
    weights : list[float]
        the allocation percentages of each strata if the allocation method is 'manual'
    mindist : float
        the minimum distance allowed between sample points
    access : SpatialVector
        a vector of LineString or MultiLineString geometries to sample near to
    layer_name : str
        the layer within 'access' to use for access buffering
    buff_inner : float
        the inner buffer around the access LineStrings, where samples should not occur
    buff_outer : float
        the outer buffer around the access LineStrings, The area in which samples must occur
    plot : bool
        whether or not to plot the output samples
    filename : str
        the output filename to write to if desired
    """
    
    if type(band) is str:
        if band not in strat_rast.bands:
            msg = "band " + str(band) + " not in given raster."
            raise ValueError(msg);

        band = strat_rast.get_band_index(band) #get band as 0-indexed integer
    else: #type(band) is int
        if band >= len(strat_rast.bands):
            msg = "0-indexed band of " + str(band) + " given, but raster only has " + str(len(raster.bands)) + " bands."
            raise ValueError(msg)

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
            raise ValueError("length of 'weights' must be the same as the number of strata.")

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
            if (mrast_band >= len(strat_rast.bands)):
                msg = "0-indexed band of " + str(band) + "given, but raster only has " + str(lend(raster.bands)) + " bands."
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

    if strat_rast.band_count != 1:
        raise ValueError("strat_raster must have a single band.")

    if not strat_rast.have_temp_dir:
        strat_rast.temp_dir = tempfile.mkdtemp()
        strat_rast.have_temp_dir = True

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
        strat_rast.temp_dir
    )

    if num_points < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, num_points))
    
    #plot new vector if requested
    if plot:
        try:
            fig, ax = plt.subplots()
            strat_rast.plot(ax, band=strat_rast.bands[0])
            title = "samples on " + strat_rast.bands[0]

            if access:
                access.plot('LineString', ax)
                title += " with access"

            ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
            ax.set_title(label=title)
            plt.show()
        except Exception as e:
            print("unable to plot output: " + str(e))

    return SpatialVector(samples)
