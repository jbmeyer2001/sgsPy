# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification by user defined quantiles
#  Author: Joseph Meyer
#  Date: September, 2025
#
# ******************************************************************************

##
# @defgroup user_quantiles quantiles
# @ingroup user_stratify

import os
import sys
import site
import tempfile
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from sgspy.utils import SpatialRaster

#ensure _sgs binary can be found
site_packages = list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0]
sys.path.append(os.path.join(site_packages, "sgspy"))
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from _sgs import quantiles_cpp, dist_cpp

GIGABYTE = 1073741824

##
# @ingroup user_quantiles
# This function conducts stratification on the raster given by generating quantile
# probabilities according to the 'quantiles' argument given by the user.
# 
# The quantiles may be defined as an integer, indicating the number of quantiles
# of equal size. Quantiles may also be defined as a list of probabilities between 0
# and 1. In the case of a raster with a single band, the quantiles may be passed directly
# to the quantiles argument as either type: int | list[float].
# 
# In the case of a multi-band raster image, the specific bands can be specified by the index
# of a list containing an equal number of quantiles as bands (list[int | list[float]). 
# 
# If not all raster bands should be stratified, specific bands can be selected in 
# the form of a dict where the key is the name of a raster band and the value is the 
# quantiles (dict[str, int | list[float]).
# 
# if the map parameter is given, an extra output band will be used which combines
# all stratifications from the bands used into an extra outpu band. A single
# value in the mapped output band corresponds to a combination a single combination
# of values from the previous bands.
# 
# The thread_count parameter specifies the number of threads which this function 
# will utilize the the case where the raster is large and may not fit in memory. If
# the full raster can fit in memory and does not need to be processed in blocks, this
# argument will be ignored. The default is 8 threads, although the optimal number
# will depend significantly on the hardware being used and may be more or less
# than 8.
# 
# the driver_options parameter is used to specify creation options for the output
# raster. See options for the Gtiff driver here: 
# https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options
# The keys in the driver_options dict must be strings, the values are converted to
# string. THe options must be valid for the driver corresponding to the filename,
# and if filename is not given they must be valid for the GTiff format, as that
# is the format used to store temporary raster files. Note that if this parameter
# is given, but filename is not and the raster fits entirely in memory, the 
# driver_options parameter will be ignored.
# 
# the eps parameter is used only if batch processing is used to calculate the quantiles
# for a raster. Quantile streaming algorithms cannot be perfectly accurate, as this
# would necessitate having the entire raster in memory at once. A good approximation
# can be made, and the error is controlled by this epsilon (eps) value.
# The Quantile streaming method is the method introduced by  Zhang et al. and utilized by MKL:
#     https://web.cs.ucla.edu/~weiwang/paper/SSDBM07_2.pdf
#     https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-summary-statistics-notes/2021-1/computing-quantiles-with-vsl-ss-method-squants-zw.html
#
# Additionally, the 'plot' parameter determines whether a histogram plot will be made of the
# stratified bands. The histogram will be the distribution of each of the raster bands which
# had stratified bands made from them, with indicators on the break values. The 'histogram_bins'
# parameter indicates the number of bands in each of the plotted histograms.
#
# The 'info' parameter, when true, prints the quantile values of each raster band.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster('rast.tif') @n
# srast = sgspy.stratify.quantiles(rast, quantiles=5)
# 
# rast = sgspy.SpatialRaster('rast.tif') @n
# srast = sgspy.stratify.quantiles(rast, quantiles=[.1, .2, .3, .5, .7], filename="srast.tif")
# 
# rast = sgspy.SpatialRaster('multi_band_rast.tif') @n
# srast = sgspy.stratify.quantiles(rast, quantiles=[5, 5, [.5, .75]], map=True)
# 
# rast = sgspy.SpatialRaster('multi_band_rast.tif') @n
# srast = sgspy.stratify.quantiles(rast, quantiles={'zq90': 5})
# 
# Parameters
# --------------------
# rast : SpatialRaster @n
#     raster data structure containing the raster to stratify @n @n
# quantiles : int | list[float] | list[int|list[float]] | dict[str,int|list[float]] @n
#     specification of the quantiles to stratify @n @n
# map : bool @n
#     whether to map the stratifiction of multiple raster bands onto a single band @n @n
# filename : str @n
#     filename to write to or ''  if no file should be written @n @n
# thread_count : int @n
#     the number of threads to use when multithreading large images @n @n
# driver_options : dict[] @n
#     the creation options as defined by GDAL which will be passed when creating output files @n @n
# eps : float @n
#     the epsilon value, controlling the error of stream-processed quantiles @n @n
# plot : optional[bool] @n
#     whether or not to plot a histogram of the values in the bands with indicators on the breaks @n @n
# histogram_bins : Optional[int] @n
#     The number of bins in the plotted histogram. @n @n
# info : Optional[bool] @n
#     when true, plot quantile values of each band after calculation @n @n
# 
# Returns
# --------------------
# a SpatialRaster object containing stratified raster bands.
def quantiles(
    rast: SpatialRaster,
    quantiles: int | list[float] | list[int|list[float]] | dict[str,int|list[float]],
    map: bool = False,
    filename: str = '',
    thread_count: int = 8,
    driver_options: dict = None,
    eps: float = .001,
    plot: Optional[bool] = None,
    histogram_bins: Optional[int] = None,
    info: Optional[bool] = None):
    
    MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow
    
    if type(rast) is not SpatialRaster:
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster")

    if type(quantiles) not in [int, list, dict]:
        raise TypeError("'quantiles' parameter must be of type int, list, or dict.")

    if type(map) is not bool:
        raise TypeError("'map' parameter must be of type bool.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if type(thread_count) is not int:
        raise TypeError("'thread_count' parameter must be of type int.")

    if type(eps) is not float:
        raise TypeError("'eps' parameter must be of type float.")

    if rast.closed:
            raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    if type(quantiles) is list and len(quantiles) < 1:
        raise ValueError("quantiles list must contain at least one element")

    if plot is not None and type(plot) is not bool:
        raise TypeError("'plot' parameter, if given, must be of type bool.")

    if histogram_bins is not None and type(histogram_bins) is not int:
        raise TypeError("'histogram_bins' parameter, if given, must be of type bool.")

    if info is not None and type(info) is not bool:
        raise TypeError("'info' parameter, if given, must be of type bool.")

    probabilities_dict = {}
    if type(quantiles) is int:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("quantiles int is for a single rast band, but the raster has {}".format(rast.band_count))

        #add quantiles to probabilities_dict
        inc = 1 / quantiles
        probabilities_dict[0] = np.array(range(1, quantiles)) / quantiles

    elif type(quantiles) is list and type(quantiles[0]) is float:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("quantiles list[float] type is for a single raster band, but the raster has {}".format(rast.band_count))

        #error check list values
        if min(quantiles) < 0:
            raise ValueError("list[float] must not contain a value less than 0")
        elif max(quantiles) > 1:
            raise ValueError("list[float] must not contain a value greater than 1")

        #add quantiles to probabilities_dict and ensure 1 and 0 are removed
        probabilities_dict[0] = quantiles
        if 0.0 in probabilities_dict[0]:
            probabilities_dict[0].remove(0.0)
        if 1.0 in probabilities_dict[0]:
            probabilities_dict[0].remove(1.0)

    elif type(quantiles) is list:
        #error checking number of raster bands
        if (len(quantiles)) != rast.band_count:
            raise ValueError("number of lists in quantiles must be equal to the number of raster bands.")

        #for each given quantiles, add it to probabilities_dict depending on type
        for i in range(len(quantiles)):
            if type(quantiles[i]) is int:
                inc = 1 / quantiles[i]
                probabilities_dict[i] = np.array(range(1, quantiles[i])) / quantiles[i]
            else: #list of float
                #for lists, error check max and min values
                if min(quantiles[i]) < 0:
                    raise ValueError("list[float] must not contain value less than 0")
                elif max(quantiles[i]) > 1:
                    raise ValueError("list[float] must not contain value greater than 1")
                probabilities_dict[i] = quantiles[i]
                if 0.0 in probabilities_dict[i]:
                    probabilities_dict[i].remove(0.0)
                if 1.0 in probabilities_dict[i]:
                    probabilities_dict[i].remove(1.0)

    else: #type dict
        for key, val in quantiles.items():
            if key not in rast.bands:
                raise ValueError("probabilities dict key must be valid band name (see SpatialRaster.bands for list of names)")
            else:
                band_num = rast.band_name_dict[key]
                if type(val) is int:
                    inc = 1 / val
                    probabilities_dict[band_num] = np.array(range(1, val)) / val
                else: #list of float
                    #for lists, error check max and min values
                    if min(val) < 0:
                        raise ValueError("list[float] must not contain value less than 0")
                    elif max(val) > 1:
                        raise ValueError("list[float] must not contain value greater than 1")
                    probabilities_dict[band_num] = val
                    if 0.0 in probabilities_dict[band_num]:
                        probabilities_dict[band_num].remove(0.0)
                    if 1.0 in probabilities_dict[band_num]:
                        probabilities_dict[band_num].remove(1.0)

    #error check max value for potential overflow error
    max_mapped_strata = int(map)
    for _, val in probabilities_dict.items():
        strata_count = len(val) + 1
        if strata_count > MAX_STRATA_VAL:
            raise ValueError("one of the quantiles given will cause an integer overflow error because the max strata number is too large.")
        max_mapped_strata = max_mapped_strata * strata_count

    if max_mapped_strata > MAX_STRATA_VAL:
        raise ValueError("the mapped strata will cause an overflow error because the max strata number is too large.")

    if thread_count < 1:
        raise ValueError("number of threads can't be less than 1.")

    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise ValueError("the key for all key/value pairs in the driver_options dict must be a string.")
            driver_options_str[key] = str(val)

    large_raster = False
    raster_size_bytes = 0
    height = rast.height
    width = rast.width
    for key, _ in probabilities_dict.items():
        pixel_size = rast.cpp_raster.get_raster_band_type_size(key)
        band_size = height * width * pixel_size
        raster_size_bytes += band_size
        if band_size >= GIGABYTE:
            large_raster = True
            break

    #if large_raster is true, the C++ function will process the raster in blocks
    large_raster = large_raster or (raster_size_bytes > GIGABYTE * 4)

    #make a temp directory which will be deleted if there is any problem when calling the cpp function
    temp_dir = tempfile.mkdtemp()
    rast.have_temp_dir = True
    rast.temp_dir = temp_dir

    #call stratify quantiles function
    [srast, quantile_vals] = quantiles_cpp(
        rast.cpp_raster, 
        probabilities_dict, 
        map, 
        filename,
        temp_dir,
        large_raster,
        thread_count,
        driver_options_str,
        eps
    )

    srast = SpatialRaster(srast)

    #now that it's created, give the cpp raster object ownership of the temporary directory
    rast.have_temp_dir = False
    srast.cpp_raster.set_temp_dir(temp_dir)
    srast.temp_dataset = filename == "" and large_raster
    srast.filename = filename

    if info:
        for band, vals in quantile_vals.items():
            print("band " + str(band) + " quantile values:")
            print(vals)
            print()
    
    if plot: 
        cpp_vector = None
        layer = ""
        bin_count = histogram_bins if histogram_bins is not None else 50

        for band, vals in quantile_vals.items():
            result = dist_cpp(rast.cpp_raster, rast.bands.index(band), cpp_vector, layer, bin_count, thread_count)
            [bins, counts] = result["population"]
            freq = counts / np.sum(counts)
            bin_size = bins[1] - bins[0]

            plt.bar(bins[0:bin_count], freq, alpha=0.5, width=bin_size, label="frequencies")
            for val in vals: plt.axvline(x=val, color='r')
            plt.legend(loc='upper right')
            plt.title(band)
            plt.show()

    metadata_info = []
    mapped_band_metadata = []
    mapped_strata_count = 1
    for band, vals in quantile_vals.items():
        name = rast.bands[band]
        strata_count = len(vals) + 1

        metadata = ["{} < {}".format(name, vals[0])]
        for i in range(1, len(vals) - 1):
            metadata.append("{} <= {} < {}".(vals[i - 1], name, vals[i]))
        metadata.append("{} <= {}".format(vals[-1], name))
        metadata_info.append(StratRasterBandMetadata(mapped=False, strata_count=strata_count, band_metadata = metadata))
        
        if map:
            mapped_band_metadata.append((srast.bands[band], strata_count))
            mapped_strata_count = mapped_strata_count * strata_count

    if map:
        metadata_info.append(StratRasterBandMetadata(mapped=True, strata_count=mapped_strata_count, mapped_band_metadata=mapped_band_metadata)

    srast.srast_metadata_info = metadata_info
    srast.is_strat_rast = True
    return srast
