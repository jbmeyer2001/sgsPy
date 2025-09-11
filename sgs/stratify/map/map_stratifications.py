# ******************************************************************************
#
#  Project: sgs
#  Purpose: map mulitiple stratification rasters
#  Author: Joseph Meyer
#  Date: September, 2025
#
# ******************************************************************************

import tempfile
from sgs.utils import SpatialRaster
from map_stratifications import map_stratifications_cpp

GIGABYTE = 1073741824
MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

def map(*args: tuple[SpatialRaster, int|str|list[int]|list[str], int|list[int]],
        filename: str = '',
        thread_count: int = 8,
        driver_options: dict = None):
    """
    this function conducts mapping on existing stratifications.

    The pre-existing stratifications are passed in the form of a raster, band, and num_stratum.
    The bands argument specifies which bands within the raster should be used, the num_stratum
    argument specifies the number of stratum within one particular band.

    the arguments are passed in the form of a tuple, of which there can be any number.
    For example, both of the following are valid:
     - map((rast1, bands1, num_stratum1))
     - map((rast1, bands1, num_stratum1), (rast1, bands2, num_stratum1))

    the raster within the tuple MUST be of type sgs.utils.SpatialRaster. 
    The bands argument MUST be: 
     - an int, specifying a single band.
     - a str, specifying a single band.
     - a list of ints, specifying the indexes of bands.
     - a list of strings, specifying the names of bands.
    
    The num_stratum argument MUST be
     - an int, if bands is an int or string, specifiying the exact number of stratum in the 
            selected band.
     - a list of ints of the same length of bands, specifying the exact number of stratum in 
            each of the indexes specified by the bands list.

    Parameters
    --------------------
    *args : tuple[SpatialRaster, int|list[int]|list[str], int|list[int]]
        tuples specifying raster bands and their number of stratifications
    filename : str
        filename to write to or '' if not file should be written
    Raises
    --------------------
    TypeError
        if one of bands or num_stratum is a list but the other is not
    ValueError
        if bands and num_stratum are both lists, but have different lengths
    ValueError
        if bands/num_stratum lists have more elements than the raster
    ValueError
        if a string within the bands argument does not exist in the raster
    ValueError
        if an int within the bands argument does not exist in the raster
    ValueError
        if the height or width of all rasters doesn't match
    ValueError
        if the maximum strata value would result in an integer overflow error
    RuntimeError (C++)
        if raster pixel type is not GDT_Float32
    RuntimeError (C++)
        if the number of mappings is large enought that it would cause integer overflow
    """

    raster_list = []
    band_lists = []
    strata_lists = []

    height = args[0][0].height
    width = args[0][0].width

    raster_size_bytes = 0
    large_raster = False
    for (raster, bands, num_stratum) in args: 
        if raster.height != height:
            raise ValueError("height is not the same across all rasters.")

        if raster.width != width:
            raise ValueError("width is not the same across all rasters.")

        #error checking on bands and num_stratum lists
        if type(bands) is list and type(num_stratum) is list and len(bands) != len(num_stratum):
            raise ValueError("if bands and num_stratum arguments are lists, they must have the same length.")
        
        if (type(bands) is list) ^ (type(num_stratum) is list): #XOR
            raise TypeError("if one of bands and num_stratum is list, the other one must be a list of the same length.")

        if type(bands) is list and len(bands) > raster.band_count:
            raise ValueError("bands list cannot have more bands than raster contains.")
            
        #helper function which checks int/str value and returns int band index
        def get_band_int(band: int|str) -> int:
            #if an int is passed, check and return
            if type(band) is int:
                if band not in range(raster.band_count):
                    raise ValueError("band {} is out of range.".format(band))
                return band

            #if a string is passed, check and return corresponding int
            if band not in raster.bands:
                msg = "band {} is not a band within the raster.".format(band)
                raise ValueError(msg)
            return raster.band_name_dict[band]

        #error checking on band int/string values
        band_list = []
        stratum_list = []
        if type(bands) is list:
            for i in range(len(bands)):
                band_int = get_band_int(bands[i])
                band_list.append(band_int)
                stratum_list.append(num_stratum[i])
                
                #check for large raster
                pixel_size = raster.cpp_raster.get_raster_band_type_size(band_int)
                band_size = height * width * pixel_size
                raster_size_bytes += band_size
                if band_size > GIGABYTE:
                    large_raster = True
        else:
            band_int = get_band_int(bands)
            band_list.append(band_int)
            stratum_list.append(num_stratum)
            
            #check for large raster
            pixel_size = raster.cpp_raster.get_raster_band_type_size(band_int)
            band_size = height * width * pixel_size
            raster_size_bytes += band_size
            if band_size > GIGABYTE:
                large_raster == true
        
        #prepare cpp function arguments
        raster_list.append(raster.cpp_raster)
        band_lists.append(band_list)
        strata_lists.append(stratum_list)

    #if any 1 band is larger than a gigabyte, or all bands together are larger than 4
    #large_raster is defined to let the C++ function know to process in blocks rather
    #than putting the entire raster into memory.
    large_raster = large_raster or (raster_size_bytes > GIGABYTE * 4)

    #error check max value for potential overflow error 
    max_mapped_strata = 1
    for strata_list in strata_lists:
        for strata_count in strata_list:
            max_mapped_strata = max_mapped_strata * strata_count
    if max_mapped_strata > MAX_STRATA_VAL:
        raise ValueError("the mapped strata will cause an overflow error because the max strata number is too large.")

    #emsire driver options keys are strings, and convert driver options vals to strings
    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise ValueError("the key for all key/value pairs in teh driver_options dict must be a string")
            driver_options_str[key] = str(val)

    temp_dir = tempfile.mkdtemp()

    #call cpp map function
    srast = SpatialRaster(map_stratifications_cpp(
        raster_list, 
        band_lists, 
        strata_lists, 
        filename, 
        large_raster,
        thread_count,
        temp_dir,
        driver_options_str
    ))

    #give srast ownership of its own temp directory
    srast.have_temp_dir = True
    srast.temp_dir = temp_dir

    return srast
