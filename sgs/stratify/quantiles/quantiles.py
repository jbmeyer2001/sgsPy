# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification by user defined quantiles
#  Author: Joseph Meyer
#  Date: September, 2025
#
# ******************************************************************************

import tempfile
import numpy as np
from sgs.utils import SpatialRaster
from quantiles import quantiles_cpp

GIGABYTE = 1073741824
MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

def quantiles(
    rast: SpatialRaster,
    num_strata: int | list[float] | list[int|list[float]] | dict[str,int|list[float]],
    map: bool = False,
    filename: str = '',
    thread_count: int = 8,
    driver_options: dict = None,
    eps: float = .001):
    """
    This function conducts stratification on the raster given by generating quantile
    probabilities according to the 'num_strata' argument given by the user.

    The quantiles may be defined as an integer, indicating the number of quantiles
    of equal size. Quantiles may also be defined as a list of probabilities between 0
    and 1. In the case of a raster with a single band, the quantiles may be passed directly
    to the num_strata argument as either type: int | list[float].

    In the case of a multi-band raster image, the specific bands can be specified by the index
    of a list containing an equal number of quantiles as bands (list[int | list[float]). 

    If not all raster bands should be stratified, specific bands can be selected in 
    the form of a dict where the key is the name of a raster band and the value is the 
    quantiles (dict[str, int | list[float]).

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to stratify
    num_strata : int | list[float] | list[int|list[float]] | dict[str,int|list[float]]
        specification of the quantiles to stratify
    map : bool
        whether to map the stratifiction of multiple raster bands onto a single band
    filename : str
        filename to write to or ''  if no file should be written

    Raises
    --------------------
    ValueError
        if number of bands required by the format/values of num_strata is inequal to the number of raster bands
    ValueError
        if quantiles specified by a list of floats (list[float]) have a maximum over 1
    ValueError
        if quantiles specified by a list of floats (list[float]) have a minimum under 0
    ValueError
        if maximum strata value would result in an integer overflow error
    RuntimeError (C++)
        if the data type of the raster is not accepted
    RuntimeError (C++)
        if num_strata is too large to the point where it would result in integer overflow
    RuntimeError (C++)
        if the number of strata required by mapping is too large to the point where it would result in integer overflow    
    """

    if type(num_strata) is list and len(num_strata) < 1:
        raise ValueError("num_strata list must contain at least one element")

    probabilities_dict = {}
    if type(num_strata) is int:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("num_strata int is for a single rast band, but the raster has {}".format(rast.band_count))

        #add quantiles to probabilities_dict
        inc = 1 / num_strata
        probabilities_dict[0] = np.array(range(1, num_strata)) / num_strata

    elif type(num_strata) is list and type(num_strata[0]) is float:
        #error check number of raster bands
        if rast.band_count != 1:
            raise ValueError("num_strata list[float] type is for a single raster band, but the raster has {}".format(rast.band_count))

        #error check list values
        if min(num_strata) < 0:
            raise ValueError("list[float] must not contain a value less than 0")
        elif max(num_strata) > 1:
            raise ValueError("list[float] must not contain a value greater than 1")

        #add quantiles to probabilities_dict and ensure 1 and 0 are removed
        probabilities_dict[0] = num_strata
        if 0.0 in probabilities_dict[0]:
            probabilities_dict[0].remove(0.0)
        if 1.0 in probabilities_dict[0]:
            probabilities_dict[0].remove(1.0)

    elif type(num_strata) is list:
        #error checking number of raster bands
        if (len(num_strata)) != rast.band_count:
            raise ValueError("number of lists in num_strata must be equal to the number of raster bands.")

        #for each given num_strata, add it to probabilities_dict depending on type
        for i in range(len(num_strata)):
            if type(num_strata[i]) is int:
                inc = 1 / num_strata[i]
                probabilities_dict[i] = np.array(range(1, num_strata[i])) / num_strata[i]
            else: #list of float
                #for lists, error check max and min values
                if min(num_strata[i]) < 0:
                    raise ValueError("list[float] must not contain value less than 0")
                elif max(num_strata[i]) > 1:
                    raise ValueError("list[float] must not contain value greater than 1")
                probabilities_dict[i] = num_strata[i]
                if 0.0 in probabilities_dict[i]:
                    probabilities_dict[i].remove(0.0)
                if 1.0 in probabilities_dict[i]:
                    probabilities_dict[i].remove(1.0)

    else: #type dict
        for key, val in num_strata.items():
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
        print("PYTHON probabilities list:")
        print(val)
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

    temp_dir = tempfile.mkdtemp()

    #call stratify quantiles function
    srast = SpatialRaster(quantiles_cpp(
        rast.cpp_raster, 
        probabilities_dict, 
        map, 
        filename,
        temp_dir,
        large_raster,
        thread_count,
        driver_options_str,
        eps
    ))

    #give srast ownership of it's own temp directory
    srast.have_temp_dir = True
    srast.temp_dir = temp_dir

    return srast
