# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratified random sampling (srs)
#  Author: Joseph Meyer
#  Date: July, 2025
#
# ******************************************************************************

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgs.utils import(
    SpatialRaster,
    SpatialVector,
    plot,
)

from strat import (
    strat_cpp,
    strat_cpp_access,
)

def strat(
    strat_rast: SpatialRaster, #TODO add band name for strat rast
    num_samples: int,
    num_strata: int,
    wrow: int = 3,
    wcol: int = 3,
    allocation: str = "prop",
    method: str = "Queinnec",
    weights: Optional[list[float]] = None,
    mindist: Optional[float] = None,
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
        the raster to sample, MUST have pixels of type 'float32' and be zero indexed
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

    Raises
    --------------------
    ValueError
        if method is not either 'random' or 'Queinnec'
    ValueError
        if allocation is not one of 'prop', 'optim', 'equal', or 'manual'
    NotImplementedError
        if 'optim' allocation is selected
    ValueError
        if 'manual' allocation is selected and 'weights' parameter is not a list
        of length num_strata which sums to 1
    ValueError
        if either wrow or wcol are less than 1 or even
    ValueError
        if an access vector is given which has more than 1 layer and layer_name parameter is not given
    ValueError
        if layer_name is not the name of a layer in the access vector
    ValueError
        if an access vector is given, and buff_inner is either not given or is less than 0
    ValueError
        if an access vector is given, and buff_outer is zero or less
    ValueError
        if an access vector is given, and buff_inner is greater than buff_outer
    ValueError
        if mindist is less than 0
    RuntimeError (C++)
        if strat_raster pixel type is not Float32
    """
    if method not in ["random", "Queinnec"]:
        raise ValueError("method must be either 'random' or 'Queinnec'.")

    if allocation not in ["prop", "optim", "equal", "manual"]:
        raise ValueError("method must be one of 'prop', 'optim', 'equal', or 'manual'.")

    if allocation == "optim":
        raise NotImplementedError("optim has not been implemented yet.")

    if allocation == "manual":
        if weights is None:
            raise ValueError("for manual allocation, weights must be given.")

        if np.sum(weights) != 1:
            raise ValueError("weights must sum to 1.")

        if len(weights) != num_strata:
            raise ValueError("length of 'weights' must be the same as the number of strata.")

    if wrow % 2 == 0 or wrow < 1:
        raise ValueError("wrow must be odd, and greater than 0.")

    if wcol % 2 == 0 or wcol < 1:
        raise ValueError("wcol must be odd, and greater then 0.")

    if access:
        if layer_name is None:
            if len(access.layers) > 1:
                raise ValueError("if there are multiple layers in the access vector, layer_name must be defined.")
            layer_name = access.layers[0]

        if buff_inner is None:
            buff_inner = 0

        if buff_outer <= 0:
            raise ValueError("buff_outer must be larger than 0.")

        if buff_inner < 0:
            raise ValueError("buff_inner can't be less than 0.")

        if buff_outer is None:
            raise ValueError("if an access vector is given, buff_outer must be defined.")

        if buff_inner >= buff_outer:
            raise ValueError("buff_outer must be greater than buff_inner.")

    if allocation != "manual":
        weights = []

    if mindist is None:
        mindist = 0

    if mindist < 0:
        raise ValueError("mindist must be greater than or equal to 0")

    if strat_rast.band_count != 1:
        raise ValueError("strat_raster must have a single band.")

    if access:
        [sample_coordinates, samples] = strat_cpp_access(
            strat_rast.cpp_raster,
            num_samples,
            num_strata,
            wrow,
            wcol,
            allocation,
            method,
            weights,
            mindist,
            access.cpp_vector,
            layer_name,
            buff_inner,
            buff_outer,
            filename,
        )

    else:
        [sample_coordinates, samples] = strat_cpp(
            strat_rast.cpp_raster,
            num_samples,
            num_strata,
            wrow,
            wcol,
            allocation,
            method,
            weights,
            mindist,
            filename,
        )

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
