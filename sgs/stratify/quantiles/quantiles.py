
# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification by user defined quantiles
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import numpy as np

from sgs.utils import SpatialRaster

from quantiles import quantiles_cpp

def quantiles(
    rast: SpatialRaster,
    num_samples: int | list[float] | list[int|list[float]] | dict[str,int|list[float]],
    map: bool = False,
    plot: bool = False,
    filename: str = ''):
    """

    """

    if type(num_samples) is int:
        passi
        #error checking and...
        #create probabilities list not including 0 and 1

    elif type(num_samples) is list and type(num_samples[0]) is float:
        pass
        #error checking and...
        #ensure probabilities list is between 0 and 1 but doesn't contain 0 or 1

    elif type(num_samples) is list:
        pass
        #error checking and...
        #for each element generate probabilities list

    else: #dype dict
        pass
        #error checking and...
        #generate probabilities list for each band specified

    #call stratify quantiles function
    strat_raster = quantiles_cpp(rast.cpp_raster, probabilities_dict, map, filename);

    #plot distribution of stratum if requested
    if plot:
        print('plotting not implemented on strat.quantiles yet')

    return SpatialRaster(strat_raster)
