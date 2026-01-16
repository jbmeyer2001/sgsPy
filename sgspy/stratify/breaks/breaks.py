# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

##
# @defgroup user_breaks breaks
# @ingroup user_stratify

import os
import sys
import tempfile
import numpy as np
from sgspy.utils import SpatialRaster

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from _sgs import breaks_cpp

GIGABYTE = 1073741824

##
# @ingroup user_breaks
# This function conducts stratification on the raster given
# according to the user defined breaks.
# 
# The breaks may be defined as a single list of ints or floats
# in the case of a raster with a single band. Or, they may be defined
# as a list of ints or floats where the index indicates the raster band.
# Or, they may be defined as a dict where the (str) key represents
# the raster band and the value is a list of ints or floats.
# 
# if the map parameter is given, an extra output band will be used which combines
# all stratifications from the previous bands used. A single value in the mapped
# output band corresponds to a single combination of values from the previous
# bands.
# 
# the filename parameter specifies an output file name. Right now the only file format
# excepted is GTiff (.tif).
# 
# the thread_count parameter specifies the number of threads which this function will 
# utilize in the case where the raster is large and may not fit in memory. If the full
# raster can fit in memory and does not need to be processed in blocks, this argument
# will be ignored. The default is 8 threads, although the optimal number will depend significantly
# on the hardware being used and my be less or more than 8.
# 
# The driver_options parameter is used to specify creation options for a the output raster.
# See options for the Gtiff driver here: https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options
# The keys in the driver_options dict must be strings, the values are converted to string.
# The options must be valid for the driver corresponding to the filename, and if filename is not given
# they must be valid for the GTiff format, as that is the format used to store temporary raster files.
# Note that if this parameter is given, but filename is not and the raster fits entirely in memory, the
# driver_options parameter will be ignored.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster("multi_band_rast.tif") @n
# srast = sgspy.stratify.breaks(rast, breaks={"band_name1": [3, 5, 11, 18]})
# 
# rast = sgspy.SpatialRaster("single_band_rast.tif") @n
# srast = sgspy.stratify.breaks(rast, breaks=[20, 40, 60, 80], filename="breaks.tif", driver_options={"COMPRESS", "LZW"}))
# 
# rast = sgspy.SpatialRaster("multi_band_rast.tif") @n
# srast = sgspy.stratify.breaks(rast, breaks={"band_name1": [3, 5, 11, 10], "band_name2": [20, 40, 60, 80]}, map=True)
# 
# rast = sgspy.SpatialRaster("multi_band_rast.tif") @n
# srast = sgspy.stratify.breaks(rast, breaks=[[3, 5, 11, 18], [40, 60, 80], [2, 5]])
# 
# Parameters
# --------------------
# rast : SpatialRaster @n
#     raster data structure containing the raster to stratify @n @n
# breaks :  list[float | list[float]] | dict[str, list[float]], @n
#     user defined breaks to stratify @n @n
# map : bool @n
#     whether to map the stratification of multiple raster bands onto a single band @n @n
# filename : str @n
#     filename to write to or '' if no file should be written @n @n
# thread_count : int @n
#     the number of threads to use when multithreading large images @n @n
# driver_options : dict[] @n
#     the creation options as defined by GDAL which will be passed when creating output files @n @n
# 
# Returns
# --------------------
# a SpatialRaster object containing stratified raster bands.
def breaks(
    rast: SpatialRaster,
    breaks: list[float | list[float]] | dict[str, list[float]],
    map: bool = False,
    filename: str = '',
    thread_count: int = 8,
    driver_options: dict = None
    ):

    MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow
    
    if type(rast) is not SpatialRaster:
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster")

    if type(breaks) not in [list, dict]:
        raise TypeError("'breaks' parameter must be of type list or dict.")

    if type(map) is not bool:
        raise TypeError("'map' parameter must be of type bool.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if type(thread_count) is not int:
        raise TypeError("'thread_count' parameter must be of type int.")

    if driver_options is not None and type(driver_options) is not dict:
        raise TypeError("'driver_options' parameter, if givne, must be of type dict.")

    if rast.closed:
            raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    breaks_dict = {}
    large_raster = False
    temp_folder = ""

    if type(breaks) is list and len(breaks) < 1:
        raise ValueError("breaks list must contain at least one element.")

    if type(breaks) is list and type(breaks[0]) is list:
        #error check number of rasters bands
        if len(breaks) != rast.band_count:
            raise ValueError("number of lists of breaks must be equal to the number of raster bands.")

        for i in range(len(breaks)):
            breaks_dict[i] = breaks[i]

    elif type(breaks) is list and type(breaks[0]) in [int, float]:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("if breaks is a single list, raster must have a single band (has {}).".format(rast.band_count))

        breaks_dict[0] = breaks

    elif type(breaks) is list:
        raise TypeError("if 'breaks' parameter is of type list, it must be filled with with values of type list, int, or float.")

    else: #breaks is a dict
        for key, val in breaks.items():
            if type(key) is not str:
                raise TypeError("if 'breaks' parameter is a dict, all keys must be of type str.")
            if type(val) is not list:
                raise TypeError("if 'breaks' parameter is a dict, all values in the key values pairs must be of type list[float].")
            if key not in rast.bands:
                raise ValueError("breaks dict key must be a valid band name (see SpatialRaster.bands for list of names)")
            
            breaks_dict[rast.band_name_dict[key]] = val

    #error check max value for potential overflow error
    max_mapped_strata = int(map)
    for _, val in breaks_dict.items():
        strata_count = len(val) + 1
        if strata_count > MAX_STRATA_VAL:
            raise ValueError("one of the breaks given will cause an integer overflow error because the max strata number is too large.")

        max_mapped_strata = max_mapped_strata * strata_count

    if max_mapped_strata > MAX_STRATA_VAL:
        raise ValueError("the mapped strata will cause an overflow error because the max strata number is too large.")    

    if thread_count < 1:
        raise ValueError("number of threads can't be less than 1.")

    #ensure driver options keys are string, and convert driver options vals to string
    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise ValueError("the key for all key/value pairs in the driver_options dict must be a string.")
            driver_options_str[key] = str(val)

    raster_size_bytes = 0
    height = rast.height
    width = rast.width
    for key, _ in breaks_dict.items():
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

    #call stratify breaks function
    srast = SpatialRaster(breaks_cpp(
        rast.cpp_raster, 
        breaks_dict, 
        map, 
        filename,
        large_raster,
        thread_count,
        temp_dir,
        driver_options_str
    ))

    #now that it's created, give the cpp raster object ownership of the temporary directory
    rast.have_temp_dir = False
    srast.cpp_raster.set_temp_dir(temp_dir)
    srast.temp_dataset = filename == "" and large_raster
    srast.filename = filename

    return srast
