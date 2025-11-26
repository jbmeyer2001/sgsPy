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

GIGABYTE = 1073741824
MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

import _sgs.map_cpp as map_cpp

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

    the filename parameter specifies an output file name. Right now the only file format
    accepted is GTiff (.tiff).

    The thread_count parameter specifies the number of threads which this function will
    utilize in the case where the raster is large an may not fit in memory. If the full
    raster can fit in memory and does not need to be processed in blocks, this argument
    will be ignored. The default is 8 threads, although the optimal number will depend
    significantly on the hardware being used and may be more or less than 8.

    the driver_options parameter is used to specifiy creation options for the output 
    raster, such as compression. See options fro GTiff driver here:
    https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options
    The keys in the driver_options dict must be strings, the values are converted to
    string. THe options must be valid for the driver corresponding to the filename,
    and if filename is not given they must be valid for the GTiff format, as that
    is the format used to store temporary raster files. Note that if this parameter
    is given, but filename is not and the raster fits entirely in memory, the 
    driver_options parameter will be ignored.

    Parameters
    --------------------
    *args : tuple[SpatialRaster, int|list[int]|list[str], int|list[int]]
        tuples specifying raster bands and their number of stratifications
    filename : str
        filename to write to or '' if not file should be written
    thread_count : int
        the number of threads to use when multithreading large images
    driver_options :
        the creation options as defined by GDAL which will be passed when creating output files

    Returns
    --------------------
    a SpatialRaster object containing a band of mapped stratifications from the input raster(s).
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
                large_raster == True
        
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

    return srast
