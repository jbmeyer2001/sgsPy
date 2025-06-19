# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import numpy as np
import matplotlib.pyplot as plt

from sgs.utils import (
    SpatialRaster,
    plot,
)

from breaks import breaks_cpp

def breaks(
    rast: SpatialRaster,
    breaks: list[int | float | list[int | float]] | dict[int, list[int|float]],
    map: bool = False,
    filename: str = ''):
    """
    This function conducts stratification on the raster given
    according to use defined breaks.

    The breaks may be defined as a single list of ints or floats
    in the case of a raster with a single band. Or, they may be defined
    as a list of ints or floats where the index indicates the raster band.
    Or, they may be defined as a dict where the (str) key represents
    the raster band and the value is a list of ints or floats.

    most of the calculation is done within the breaks_cpp function 
    which can be found in sgs/stratify/breaks/breaks.cpp/

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to stratify
    breaks : list[int | float | list[int|float]] | dict[int, list[int|float]]
        user defined breaks to stratify
    map : bool
        whether to map the stratification of multiple raster layers onto a single layer
    filename : str
        filename to write to or '' if no file should be written

    Raises
    --------------------
    ValueError
        if number of bands required by the size of the parameter 'breaks' is inequal to the number of raster bands
    ValueError
        if a break contains a value less than the minimum in the corresponding raster band
    ValueError
        if a break contains a value greater than the maximum in the corresponding raster band
    """

    #TODO breaks_dict and add dict as possible breaks parameter type
    breaks_dict = {}
    if type(breaks) is list and type(breaks[0]) is list:
        #error check number of rasters bands
        if len(breaks) != rast.band_count:
            raise ValueError("number of lists of breaks must be equal to the number of raster bands.")

        for i in range(len(breaks)):
            breaks_dict[i] = breaks[i]

    elif type(breaks) is list: #type(breaks[0]) is int or float
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("if breaks is a single list, raster must have a single band (has {}).".format(rast.band_count))

        breaks_dict[0] = breaks

    else: #breaks is a dict
        for key, val in breaks.items():
            if key in rast.bands:
                breaks_dict[rast.band_name_dict[key]] = val
            else:
                raise ValueError("breaks dict key must be a valid band name (see SpatialRaster.bands for list of names)")

    #call stratify breaks function
    strat_raster = breaks_cpp(rast.cpp_raster, breaks_dict, map, filename)

    #plot distribution of breaks if requested
    if plot:
        print('plotting not implemented on strat.breaks yet')

    return SpatialRaster(strat_raster)

