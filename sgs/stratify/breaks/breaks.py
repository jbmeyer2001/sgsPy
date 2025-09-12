# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import tempfile
import numpy as np
from sgs.utils import SpatialRaster
from breaks import breaks_cpp

GIGABYTE = 1073741824
MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

def breaks(
    rast: SpatialRaster,
    breaks: list[float | list[float]] | dict[str, list[float]],
    map: bool = False,
    filename: str = '',
    thread_count: int = 8,
    driver_options: dict = None
    ):
    """
    This function conducts stratification on the raster given
    according to the user defined breaks.

    The breaks may be defined as a single list of ints or floats
    in the case of a raster with a single band. Or, they may be defined
    as a list of ints or floats where the index indicates the raster band.
    Or, they may be defined as a dict where the (str) key represents
    the raster band and the value is a list of ints or floats.

    if the map parameter is given, an extra output band will be used which combines
    all stratifications from the previous bands used. A single value in the mapped
    output band corresponds to a single combination of values from the previous
    bands.

    the filename parameter specifies an output file name. Right now the only file format
    excepted is GTiff (.tif).

    the thread_count parameter specifies the number of threads which this function will 
    utilize in the case where the raster is large and may not fit in memory. If the full
    raster can fit in memory and does not need to be processed in blocks, this argument
    will be ignored. The default is 8 threads, although the optimal number will depend significantly
    on the hardware being used and my be less or more than 8.

    The driver_options parameter is used to specify creation options for a the output raster.
    See options for the Gtiff driver here: https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options
    The keys in the driver_options dict must be strings, the values are converted to string.
    The options must be valid for the driver corresponding to the filename, and if filename is not given
    they must be valid for the GTiff format, as that is the format used to store temporary raster files.
    Note that if this parameter is given, but filename is not and the raster fits entirely in memory, the
    driver_options parameter will be ignored.

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to stratify
    breaks :  list[float | list[float]] | dict[str, list[float]],
        user defined breaks to stratify
    map : bool
        whether to map the stratification of multiple raster bands onto a single band
    filename : str
        filename to write to or '' if no file should be written
    thread_count : int
        the number of threads to use when multithreading large images
    driver_options : dict[]
        the creation options as defined by GDAL which will be passed when creating output files
        """

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
                raise ValueError("the key for all key/value pairs in the driver_options dict must be a string")
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

    #if a VRT dataset will be used, make a temp directory to hold it's files
    temp_dir = tempfile.mkdtemp()

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

    #give srast ownership of it's own temp directory
    srast.have_temp_dir = True
    srast.temp_dir = temp_dir

    return srast
