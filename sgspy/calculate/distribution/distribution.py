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
# This function calculates the distribution of a given raster band.
# The distribution is plotted by default, it is also returned.
#
# A number of bins can be given (the default is 50), the bins will be
# sized equally depending on the minimum and maximum values within
# the raster band. If the raster is a multi-band raster, the 'band'
# parameter must be passed.
#
# Additionally, a SpatialVector may be passed to the samples parameter
# which must contain a sample plot network (only geometry types Point or 
# MultiPoint). If a SpatialVector is passed this way, the distribution
# of the samples will also be calculated, using the same bins as the layer distribution.
# if the SpatialVector has more than one layer, the 'layer' parameter must
# be passed.
#
# This function returns a dict, the keys of the dict will be either 'population' or
# 'sample'. If the 'samples' parameter is not given, only the 'population' key will exist.
# If the 'samples' parameter IS given, both 'population' and 'sample' will exist as keys.
# The values of the dict will be two arrays, one is the array of bins of type float64, the
# other is the array of counts of type int64. The array of bins will be one value longer
# than the array of counts, as it countains not just the minimum value but also the maximum value.
#
# For example, if the minimum pixel value in the raster is 1, and the maximum pixel value in the
# raster is 3, 3 and the bin size is 4, then the bin array would be [1.0, 1.5, 2.0, 2.5, 3.0]. 
# The array of counts would then be of length 4, where the first element (index 0) would be
# the count of pixels which have a pixel value that satisfies the following 1.0 <= val < 1.5.
# The second element (index 1) would then be 1.5 <= val < 2.0, and so on until the final 
# element with 2.5 <= val <= 3. This is typical, and matches the way that the numpy 
# histogram function works.
#
# The 'threads' parameter specifies the number of threads which this function will utilize.
# The default value is 8. The optimal number will depend significantly on both the
# hardware being used and the raster data, and may be more or less than 8.
#
# Examples
# --------------------
# rast = sgspy.SpatialRaster("rast.tif") #single band raster @n
# sgspy.calculate.distribution(rast) #plotting is default
#
# rast = sgspy.SpatialRaster("mraster.tif") @n
# sgspy.calculate.distribution(rast, band='pzabove2', bins=10)
#
# rast = sgspy.SpatialRaster("mraster.tif") @n
# samples = sgspy.sample.clhs(rast, num_samples=250) @n
# sgspy.calculate.distribution(rast, band='zsd', samples=samples)
#
# rast = sgspy.SpatialRaster("mraster.tif") @n
# samples = sgpsy.sample.clhs(rast, num_samples=240) @n
# result = sgspy.calculate.distribution(rast, band='zq90', samples=samples, bins=100, plot=False) @n
# [pop_bins, pop_counts] = result['population'] @n
# [samp_bins, samp_counts] = result['sample']
#
# Parameters
# --------------------
# rast : SpatialRaster @n
#   raster data structure containing input raster band @n @n
# band : Optional[str | int] @n
#   the band to use within 'rast' if 'rast' has more than one band @n @n
# samples : Optional[SpatialVector]
#   the option data structure containing an input sample network @n @n
# layer : Optional[str] @n
#   the layer name of the sample network if 'samples' contains more than one layer @n @n
# bins : int @n 
#   the number of bins in the histogram distribution @n @n 
# threads : int @n 
#   the number of threads to use @n @n 
# plot : bool @n 
#   whether to plot a histogram of the distribution @n @n 
#
# Returns
# --------------------
# a dict[str, (array, array)] containing the bin values and counts, with dict keys 'population' and potentially 'sample'
def distribution(
    rast: SpatialRaster,
    band: Optional[str | int] = None,
    samples: Optional[SpatialVector] = None,
    layer: Optional[str] = None,
    bins: int = 50,
    threads: int = 8,
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
