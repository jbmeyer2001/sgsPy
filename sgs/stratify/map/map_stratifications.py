# ******************************************************************************
#
#  Project: sgs
#  Purpose: map mulitiple stratification rasters
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

from sgs.utils import SpatialRaster

from map_stratifications import map_stratifications_cpp

def map(*args: tuple[SpatialRaster, int|str|list[int]|list[str], int|list[int]],
        filename: str = '',
        plot: bool = False):
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
     - a lsit of strings, specifying the names of bands.
    
    The num_stratum argument MUST be
     - an int, if bands is an int or string, specifiying the exact number of stratum in the 
            selected band.
     - a list of ints of the same length of bands, specifying the exact number of stratum in 
            each of the indexes specified by the bands list.

    Parameters
    --------------------
    *args : tuple[SpatialRaster, int|list[int]|list[str], int|list[int]]
        tuples specifying raster bands and their number of stratifications
    plot : bool
        whether to plot the distribution or not
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
    """
    #TODO add cpp runtime errors

    raster_list = []
    band_lists = []
    stratum_lists = []

    for (raster, bands, num_stratum) in args: 
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
                raise ValueError("bands {} is not a band within the raster.".format(band))
            return raster.band_name_dict[band]


        #error checking on band int/string values
        band_list = []
        stratum_list = []
        if type(bands) is list:
            for i in range(len(bands)):
                band_list.append(get_band_int(bands[i]))
                stratum_list.append(num_stratum[i])
        else:
            band_list.append(get_band_int(bands))
            stratum_list.append(num_stratum)
        
        #prepare cpp function arguments
        raster_list.append(raster.cpp_raster)
        band_lists.append(band_list)
        stratum_lists.append(stratum_list)

    #call cpp map function
    mapped_raster = map_stratifications_cpp(raster_list, band_lists, stratum_lists, filename)

    #plot distribtion of stratum if requested
    if plot:
        print('plotting not implemented on strat.map yet')

    return SpatialRaster(mapped_raster)
