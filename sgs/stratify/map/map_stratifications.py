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
        plot: bool = false):
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
    bands_list = []
    stratums_list = []

    for arg in args:
        raster, bands, num_stratum = args

        #error checking on bands and num_stratum lists
        if type(bands) is list and type(num_stratum) is list and len(bands) != len(num_stratum):
            raise ValueError("if bands and num_stratum arguments are lists, they must have the same length.")
        
        if type(bands) is list ^ type(num_stratum) is list: #XOR
            raise TypeError("if one of bands and num_stratum is list, the other one must be a list of the same length.")

        if type(bands) is list and len(bands) > raster.band_count:
            raise ValueError("bands list cannot have more bands than raster contains.")

        #error checking on bands string values
        if type(bands[0]) is str:
            for i in range(len(bands))
                if (bands[i] not in raster.bands):
                    raise ValueError("{} is not a band within the raster".format(bands[i]))

                #convert from string to int
                bands[i] = raster.band_name_dict[bands[i]]

        #error checking on bands int values
        if type(bands[i]) is int:
            for i in range(len(bands))
                if bands[i] is not in range(raster.band_count):
                    raise ValueError("band {} is out of range.".format(bands[i])) 

        #add raster, bands, and stratums to arguments (passed to c++ function)
        raster_list.append(raster)
        bands_list.append(bands)
        stratums_list.append(stratums)

    #call cpp map function
    mapped_raster = map_stratifications_cpp(raster_list, bands_list, stratums_list, filename)

    #plot distribtion of stratum if requested
    if plot:
        print('plotting not implemented on strat.map yet')

    return SpatialRaster(mapped_raster)
