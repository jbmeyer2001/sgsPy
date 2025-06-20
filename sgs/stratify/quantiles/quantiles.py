
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

    probabilities_dict = {}
    if type(num_samples) is int:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("num_samples int is for a single rast band, but the raster has {}".format(rast.band_count))

        #add quantiles to probabilities_dict
        inc = 1 / num_samples
        probabilities_dict[0] = np.arange(inc, 1, inc)

    elif type(num_samples) is list and type(num_samples[0]) is float:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("num_samples list[float] type is for a single raster band, but the raster has {}".format(rast.band_count))

        #error check list values
        if min(num_samples) < 0:
            raise ValueError("list[float] must not contain a value less than 0")
        elif max(num_samples) > 1:
            raise ValueError("list[float] must not contain a value greater than 1")

        #add quantiles to probabilities_dict and ensure 1 and 0 are removed
        probabilities_dict[0] = num_samples
        probabilities_dict[0].remove(0.0)
        probabilities_dict[0].remove(1.0)

    elif type(num_samples) is list:
        #error checking number of raster bands
        if (len(num_samples)) != rast.band_count:
            raise ValueError("number of lists in num_samples must be equal to the number of raster bands.")

        #for each given num_sample, add it to probabilities_dict depending on type
        for i in range(len(probabilities)):
            if type(probabilities[i]) is int:
                inc = 1 / num_samples[probabilities[i]]
                probabilities_dict[i] = np.arange(inc, 1, inc)
            else: #list of float
                #for lists, error check max and min values
                if min(probabilities[i]) < 0:
                    raise ValueError("list[float] must not contain value less than 0")
                elif max(probabilities[i]) > 1:
                    raise ValueError("list[float] must not contain value greater than 1")
                probabilities_dict[i] = probabilities[i]
                probabilities_dict[i].remove(0.0)
                probabilities_dict[i].remove(1.0)

    else: #type dict
        for key, val in num_samples.items():
            if key not in rast.bands:
                raise ValueError("probabilities dict key must be valid band name (see SpatialRaster.bands for list of names)")
            else:
                band_num = rast.band_name_dict[key]
                if type(val) is int:
                    inc = 1 / val
                    probabilities_dict[band_num] = np.arange(inc, 1, inc)
                else: #list of float
                    #for lists, error check max and min values
                    if min(val) < 0:
                        raise ValueError("list[float] must not contain value less than 0")
                    elif max(val) > 1:
                        raise ValueError("list[float] must not contain value greater than 1")
                    probabilities_dict[band_num] = val
                    probabilities_dict[band_num].remove(0.0)
                    probabilities_dict[band_num].remove(1.0)

    #call stratify quantiles function
    strat_raster = quantiles_cpp(rast.cpp_raster, probabilities_dict, map, filename);

    #plot distribution of stratum if requested
    if plot:
        print('plotting not implemented on strat.quantiles yet')

    return SpatialRaster(strat_raster)
