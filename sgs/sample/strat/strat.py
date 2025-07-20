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
    TODO add documentation 
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

        if len(weights) != num_samples:
            raise ValueError("length of 'weights' must be the same as the number of samples.")

    #is this okay?
    if wrow % 2 == 0 or wrow < 1:
        raise ValueError("wrow must be odd, and greater than 0.")

    #is this okay?
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
        [sample_coordinates, sample_points] = strat_cpp_access(
            strat_rast.cpp_raster,
            num_samples,
            num_strata,
            wrow,
            wcol,
            allocation,
            method,
            weights,
            mindist,
            access,
            layer_name,
            buff_inner,
            buff_outer,
            filename,
        )

    else:
        [sample_coordinates, sample_points] = strat_cpp(
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
    #TODO do this in try/catch so that error won't cause 
    #sampling to be thrown out
    if plot:
        fig, ax = plt.subplots()
        #TODO let user know which band is being printed
        strat_rast.plot(ax, band=strat_rast.bands[0])
        if access:
            access.plot('LineString', ax)
        ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
        plt.show()

    return sample_points
