# ******************************************************************************
#
#  Project: sgs
#  Purpose: calculate probability distribution of raster
#  Author: Joseph Meyer
#  Date: January, 2026
#
# ******************************************************************************

##
# @defgroup user_distribution distribution
# @ingroup user_calculate

import os
import sys
import site
from sgspy.utils import SpatialRaster, SpatialVector
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

#ensure _sgs binary can be found
site_packages = list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0]
sys.path.append(os.path.join(site_packages, "sgspy"))
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from _sgs import dist_cpp

##
# @ingroup user_distribution
#
# DOCUMENTATION
def distribution(
    rast: SpatialRaster,
    band: Optional[str | int] = None,
    samples: Optional[SpatialVector] = None,
    layer: Optional[str] = None,
    bins: int = 50,
    threads: int = 8,#this is automatically set to 1 on the back end right now
    plot: bool = True):

    if type(rast) is not SpatialRaster:
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster.")

    if band is not None and type(band) not in [str, int]:
        raise TypeError("'band' parameter, if given, must be of type int or string.")

    if samples is not None and type(samples) is not SpatialVector:
        raise TypeError("'samples' parameter, if given, must be of type sgspy.SpatialVector.")

    if layer is not None and type(layer) is not str:
        raise TypeError("'layer' parameter, if given, must be of type str.")

    if type(bins) is not int:
        raise TypeError("'bins' parameter must be of type int.")

    if type(plot) is not bool:
        raise TypeError("'plot' parameter must be of type bool.")

    if band is None and len(rast.bands) > 1:
        raise ValueError("'If the raster has more than 1 band, the 'band' parameter must be given.")

    if type(band) is int and band < 0:
        raise ValueError("'If the the 'band' parameter is of type int, it must be greater than or equal to 0, indicating the (zero-indexed) band index.")

    if type(band) is int and band >= len(rast.bands):
        raise ValueError("'The 'band' parameter is too large to specify one of the (zero-indexed) band indices.")

    if type(band) is str and band not in rast.bands:
        raise ValueError("'If the 'band' parameter is of type str, it must match one of the bands in the 'rast' SpatialRaster object.")
    
    band = rast.get_band_index(band)

    if samples:
        if layer is None:
            if len(samples.layers) > 1:
                raise ValueError("If there are multiple layers in the 'samples' vector, the 'layer' parameter must be passed.")

            layer = samples.layers[0]
          
        if layer not in samples.layers:
            raise ValueError("Layer specified by 'layer' parameter must be a layer in the 'samples' SpatialVector object.")

        cpp_vector = samples.cpp_vector
    else:
        cpp_vector = None
        layer = ""

    if bins < 1:
        raise ValueError("'bins' parameter must be 1 or greater.")

    #the reason why there is a cpp function written to do this, rather than just using
    #numpys histogram function is because numpys histogram function requires that all
    #the data be in a numpy array (in memory) at once, and on very large raster images 
    #this is not possible.
    result = dist_cpp(rast.cpp_raster, band, cpp_vector, layer, bins, threads)

    if plot:
        [pop_bins, pop_counts] = result["population"]
        pop_freq = pop_counts / np.sum(pop_counts)
        bin_size = pop_bins[1] - pop_bins[0]

        #essentially a histogram chart. The pyplot hist function calculates the histogram on raw data so can't
        #be used in this instance.
        plt.bar(pop_bins[0:bins], pop_freq, alpha=0.5, width=bin_size, label="population frequencies")

        if "sample" in result:
            [_, samp_counts] = result["sample"]
            samp_freq = samp_counts / np.sum(samp_counts)

            #sample bins match the population bins
            plt.bar(pop_bins[0:bins], samp_freq, alpha=0.5, width=bin_size, label="sample frequenceis")

        plt.legend(loc="upper right")
        plt.title(rast.bands[band])
        plt.show()

    return result
