# ******************************************************************************
#
#  Project: sgs
#  Purpose: map mulitiple stratification rasters
#  Author: Joseph Meyer
#  Date: September, 2025
#
# ******************************************************************************

##
# @defgroup user_map map
# @ingroup user_stratify

import os
import sys
import site
import tempfile
from sgspy.utils import SpatialRaster

#ensure _sgs binary can be found
site_packages = list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0]
sys.path.append(os.path.join(site_packages, "sgspy"))
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from _sgs import map_cpp

GIGABYTE = 1073741824

##
# @ingroup user_map
# This function conducts mapping on existing stratifications.
# 
# The pre-existing stratifications are passed in the form of a raster, band.
# The bands argument specifies which bands within the raster should be used.
# 
# **IMPORTANT**
# If the strat raster IS NOT the return value of one of the other sgspy
# stratification functions, an additional argument MUST be passed: the number of
# strata. This is because during stratification this number is automatically stored.
# The num_strata argument, if required, should be set to the value of the largest strata + 1. 
# For example if the strata are [0, 1, 2, 3, 4] then num_strata should be 5. 
# If the strata are [1, 2, 4] then num_strata should still be 5. If the strata are [0, 1, 2, 3] then
# num_strata should be 4.
# 
# the arguments are passed in the form of a tuple, and there can be any number of tuples passed.
# For example, the following are valid:
#  - map((srast1, bands1))
#  - map((srast1, bands1), (rast1, bands2))
#  - map((srast1, bands1, num_strata1)) #if rast1 is NOT the return value of an sgspy stratification function.
#  - map((srast1, bands1, num_strata1), (srast2, bands2, num_strata2)) #if rast1 and rast2 are NOT the return values of sgspy stratification functions.
# 
# the raster within the tuple MUST be of type sgs.utils.SpatialRaster. 
# The bands argument MUST be: 
#  - an int, specifying a single band in the strat raster.
#  - a str, specifying a single band in the strat raster.
#  - a list of ints or strs, each specifying a band in the strat raster.
# 
# The num_strata argument, if required, MUST be:
#  - an int, if bands argument is an int or string, specifiying the max strata value in a specific band.
#  - a list of ints, if bands argument is a list, each entry specifying the max strata value in that band.
# 
# the filename parameter specifies an output file name. Right now the only file format
# accepted is GTiff (.tiff).
# 
# The thread_count parameter specifies the number of threads which this function will
# utilize in the case where the raster is large an may not fit in memory. If the full
# raster can fit in memory and does not need to be processed in blocks, this argument
# will be ignored. The default is 8 threads, although the optimal number will depend
# significantly on the hardware being used and may be more or less than 8.
# 
# the driver_options parameter is used to specifiy creation options for the output 
# raster, such as compression. See options fro GTiff driver here:
# https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options
# The keys in the driver_options dict must be strings, the values are converted to
# string. THe options must be valid for the driver corresponding to the filename,
# and if filename is not given they must be valid for the GTiff format, as that
# is the format used to store temporary raster files. Note that if this parameter
# is given, but filename is not and the raster fits entirely in memory, the 
# driver_options parameter will be ignored.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster("rast.tif") @n
# breaks = sgspy.stratify.breaks(rast, breaks={'zq90': [3, 5, 11, 18], 'pzabove2]: [20, 40, 60, 80]}) @n
# quantiles = sgspy.stratify.quantiles(rast, quantiles={'zsd': 25}) @n
# srast = sgspy.stratify.map((breaks, ['strat_zq90', 'strat_pzabove2']), (quantiles, 'strat_zsd'))
# 
# rast = sgspy.SpatialRaster("rast.tif") @n
# inventory = sgspy.SpatialVector("inventory_polygons.shp") @n
# breaks = sgspy.stratify.breaks(rast, breaks={'zq90': [3, 5, 11, 18], 'pzabove2]: [20, 40, 60, 80]}) @n
# poly = sgspy.stratify.poly(rast, inventory, attribute="NUTRIENTS", layer_name="inventory_polygons", features=['poor', 'medium', 'rich']) @n
# srast = sgspy.stratify.map((breaks, [0, 1]), (poly, 0), filename="mapped_srast.tif", driver_options={"COMPRESS", "LZW"})
#
# #Some pre-existing strat raster(s) @n
# breaks = sgspy.SpatialRaster("breaks.tif") @n
# quantiles = sgspy.SpatialRaster("quantiles.tif") @n
# #give an extra argument, the number of strata, for each strat raster, because sgspy does not know from previously creating that sraster. @n
# mapped = sgspy.stratify.map((breaks, ['strat_zq90', 'strat_pzabove2'], [5, 5]), (quantiles, 'strat_zsd', 25))
#
# #Another example with pre-existing strat raster(s)
# rast = sgspy.SpatialRaster("rast.tif") @n
# breaks = sgspy.SpatialRaster("breaks.tif") @n
# poly = sgspy.SpatialRaster("poly.tif") @n
# #give an extra argument, the number of strata, for each strat raster, because sgspy does not know from previously creating that sraster. @n
# srast = sgspy.stratify.map((breaks, [0, 1], [5, 5]), (poly, 0, 3))
#
# Parameters
# --------------------
# *args : tuple[SpatialRaster, int|list[int]|list[str], Optional[int|list[int]]] @n
#     tuples specifying raster bands and their number of stratifications @n @n
# filename : str @n
#     filename to write to or '' if not file should be written @n @n
# thread_count : int @n
#     the number of threads to use when multithreading large images @n @n
# driver_options : dict[str]  @n
#     the creation options as defined by GDAL which will be passed when creating output files @n @n
# 
# Returns
# --------------------
# a SpatialRaster object containing a band of mapped stratifications from the input raster(s).
def map(*args: tuple[SpatialRaster, int|str|list[int]|list[str], Optional[int|list[int]]],
        filename: str = '',
        thread_count: int = 8,
        driver_options: dict = None):
            
    MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if type(thread_count) is not int:
        raise TypeError("'thread_count' parameter must be of type int.")

    if driver_options is not None and type(driver_options) is not dict:
        raise TypeError("'driver_options' parameter, if given, must be of type dict.")

    raster_list = []
    band_lists = []
    strata_lists = []

    height = None
    width = None

    raster_size_bytes = 0
    large_raster = False
    mapped_band_metadata = []
    mapped_strata_count = 1
    for (raster, bands, num_strata) in args:
        #error checking on input args
        if type(raster) is not SpatialRaster:
            raise TypeError("first value in each tuple argument, which represents the strat raster, must be of type sgspy.SpatialRaster.")

        if type(bands) not in [int, str, list]:
            raise TypeError("second value in each tuple argument, which represents the band(s) to use within the strat raster, must be of type int, str, or list.")

        if !raster.is_strat_rast and num_strata is None:
            raise TypeError("if one of the strat rasters specified is not the return value of running one of the sgspy stratification functions, the additional 'num_strata' argument must be given.")

        if raster.is_strat_rast:
            num_strata = None

        if num_strata is not None and type(num_strata) not in [int, list]:
            raise TypeError("if the num_strata argument is required, it must be of type int or list.")

        if raster.closed:
            raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

        if not height:
            height = raster.height
        if raster.height != height:
            raise ValueError("height is not the same across all rasters.")

        if not width:
            width = raster.width
        if raster.width != width:
            raise ValueError("width is not the same across all rasters.")

        if num_strata is not None and type(bands) is list and type(num_strata) is list and len(bands) != len(num_strata):
            raise ValueError("if bands and num_strata arguments are lists, they must have the same length.")
        
        if num_strata is not None and ((type(bands) is list) ^ (type(num_strata) is list)): #XOR
            raise TypeError("if the 'num_strata' argument is given, and one of the 'bands' argument and 'num_strata' argument is a list, the other one must also be a list of the same length.")

        if type(bands) is list and len(bands) > raster.band_count:
            raise ValueError("bands list cannot have more bands than raster contains.")
            
        #helper function which checks int/str value and returns int band index
        def get_band_int(band: int|str) -> int:
            #if an int is passed, check and return
            if type(band) is int:
                if band not in range(raster.band_count):
                    raise ValueError("band {} is out of range.".format(band))

            #if a string is passed, check and return corresponding int
            elif type(band) is str:
                if band not in raster.bands:
                    msg = "band {} is not a band within the raster.".format(band)
                raise ValueError(msg)

            else:
                raise TypeError("if the band argument is a list, every value within it must be of type int or str.")
            return raster.get_band_index(band)

        #add raster, bands, and num_strata to lists which will be passed to the C++ function
        if type(bands) is list:
            for i in range(len(bands)):
                band_int = get_band_int(bands[i])
                band_list.append(band_int)
                if raster.is_strat_rast:
                    num_strata = raster.srast_metadata_info[band_int].get_num_strata()
                    strata_count_list.append(num_strata)
                else:
                    strata_count_list.append(num_strata[i])
                
                #check for large raster
                pixel_size = raster.cpp_raster.get_raster_band_type_size(band_int)
                band_size = height * width * pixel_size
                raster_size_bytes += band_size
                if band_size > GIGABYTE:
                    large_raster = True
        else:
            band_int = get_band_int(bands)
            band_list.append(band_int)
            if raster.is_strat_rast:
                num_strata = raster.srast_metadata_info[band_int].get_num_strata()
                strata_count_list.append(num_strata)
            else:
                strata_count_list.append(num_strata)
            
            #check for large raster
            pixel_size = raster.cpp_raster.get_raster_band_type_size(band_int)
            band_size = height * width * pixel_size
            raster_size_bytes += band_size
            if band_size > GIGABYTE:
                large_raster == True
        
        #prepare cpp function arguments
        raster_list.append(raster.cpp_raster)
        band_lists.append(band_list)
        strata_lists.append(strata_count_list)

    #if any 1 band is larger than a gigabyte, or all bands together are larger than 4, large_raster is true
    #
    #large_raster is defined to let the C++ function know to process using an in-memory dataset (MEM) or 
    #another virtual type (VRT)
    large_raster = large_raster or (raster_size_bytes > GIGABYTE * 4)

    #error check max value for potential overflow error 
    max_mapped_strata = 1
    for strata_list in strata_lists:
        for strata_count in strata_list:
            max_mapped_strata = max_mapped_strata * strata_count
    if max_mapped_strata > MAX_STRATA_VAL:
        raise ValueError("the mapped strata will cause an overflow error because the max strata number is too large.")

    #ensire driver options keys are strings, and convert driver options vals to strings
    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise ValueError("the key for all key/value pairs in teh driver_options dict must be a string")
            driver_options_str[key] = str(val)

    #make a temp directory which will be deleted if there is any problem when calling the cpp function
    temp_dir = tempfile.mkdtemp()
    args[0][0].have_temp_dir = True
    args[0][0].temp_dir = temp_dir

    #call cpp map function
    srast = SpatialRaster(map_cpp(
        raster_list, 
        band_lists, 
        strata_lists, 
        filename, 
        large_raster,
        thread_count,
        temp_dir,
        driver_options_str
    ))

    #now that it's created, give the cpp raster object ownership of the temporary directory
    args[0][0].have_temp_dir = False
    srast.cpp_raster.set_temp_dir(temp_dir)
    srast.temp_dataset = filename == "" and large_raster
    srast.filename = filename

    mapped_band_metadata = []
    mapped_strata_count = 1
    for i in len(range(raster_list)):
        rast = raster_list[i]
        bands = band_lists[i]
        stratas = strata_lists[i]

        for j in len(range(bands))
            band = bands[j]
            strata = stratas[j]

            mapped_band_metadata.append(rast.bands[band], strata)
            mapped_strata_count = mapped_strata_count * strata

    srast.srast_metadata_info = [
        StratRasterBandMetadata(mapped=True, strata_count=mapped_strata_count, mapped_band_metadata=mapped_band_metadata)
    ]
    srast.is_strat_rast = True
    
    return srast
