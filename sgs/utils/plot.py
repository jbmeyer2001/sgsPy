import numpy as np
import matplotlib.pyplot as plt

def arrange_bands_from_list(raster, bands_list):
    """
    Used by plot_raster function. Converts all bands in initial list to int indexes
    using SpatialRaster.get_band_index().

    Parameters
    --------------------
    bands_list: list
        a list containing str or int variables specify bands

    Raises
    --------------------
    ValueError:
        if the number of bands in the band list is not 1 or 3
    """
    num_items = len(bands_list)
    if num_items != 1 and num_items != 3:
        raise ValueError("number of bands must be either 1 (for scalar iamges) or 3 (for RGB images).")

    for i in range(num_items):
        bands_list[i] = raster.get_band_index(bands_list[i])

    return bands_list

def arrange_bands_from_dict(raster, bands_dict):
    """
    Used by plot_raster function. Converts a dict which specifies RGB values
    into a list of indexes using SpatialRaster.get_band_index().

    Parameters
    --------------------
    bands_dict: dict
        dict specifying 'red', 'green', and 'blue' bands

    Raises
    --------------------
    ValueError:
        if the bands dict does not have three items, or if the items have incorrect names
    """
    if len(bands_dict) != 3 or 'red' not in bands_dict or 'green' not in bands_dict or 'blue' not in bands_dict:
        raise ValueError("if bands is a dict, it must to have three items with the keys 'red', 'green', and 'blue'.")

    return [
        raster.get_band_index(bands_dict["red"]),
        raster.get_band_index(bands_dict["green"]),
        raster.get_band_index(bands_dict["blue"])
    ]

def plot_raster(raster, ax, target_width, target_height, bands=None, **kwargs):
    """
    Plots the specified bands using matplotlib.pyplot.imshow function.

    Parameters
    --------------------
    raster : SpatialRaster
        raster to plot
    ax : matplotlib axis
        the axis to plot the image on
    target_width : int
        maximum width in pixels for the image (after downsampling)
    target_height : int
        maximum height in pxeils for the image (after downsampling)
    bands (optional) : int or str or list or dict
        specification of which bands to plot
    **kwargs (optional)
        any parameters which may be passed to matplotlib.pyplot.imshow

    Raises
    --------------------
    TypeError:
        if 'bands' is not of type int, str, list, or dict
    """
    #get bands argument as list of int
    if bands is None:
        bands = arrange_bands_from_list(raster, [*range(raster.band_count)])
    elif type(bands) == list:
        bands = arrange_bands_from_list(raster, bands)
    elif type(bands) == dict:
        bands = arrange_bands_from_dict(raster, bands)
    elif type(bands) in [str, int]: 
        bands = [raster.get_band_index(bands)]
    else:
        raise TypeError("'bands' parameter must be of type None, list, dict, str, or int.")

    #calculate downsampled resolution and get downsampled raster
    #for info on downsample resolution calculation:
    #https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec
    target_downscaling_factor = min(raster.width / target_width, raster.height / target_height)
    if (target_downscaling_factor <= 2 / 1.2):
        downsampled_width = raster.width
        downsampled_height = raster.height
    elif (target_downscaling_factor <= 4 / 1.2):
        downsampled_width = int(raster.width / 2)
        downsampled_height = int(raster.height / 2)
    elif (target_downscaling_factor <= 8 / 1.2):
        downsampled_width = int(raster.width / 4)
        downsampled_height = int(raster.height / 4)
    else:
        downsamlped_width = int(raster.width / 8)
        downsampled_heigth = int(raster.height / 8)
    arr = np.asarray(
        raster.cpp_raster.get_raster_as_memoryview(downsampled_width, downsampled_height),
        copy=False
    )

    #get raster origin and raster extent
    extent = (raster.xmin, raster.xmax, raster.ymin, raster.ymax) #(left, right, top, bottom)

    #add image to matplotlib
    display_arr = np.moveaxis(arr[bands, :, :], 0, 2)
    ax.imshow(display_arr, origin='upper', extent=extent, **kwargs)

def plot_vector(vector, ax, geomtype, layer, **kwargs):
    """
    Plots the specified layer using matplotlib.pyplot.plot.
    The parameter give by geomtype must be one of:
    'Point', 'MultiPoint', 'LineString', 'MultiLineString'.

    The layer must contain only geometries of type Point and
    MultiPoint in the case where 'Point' or 'MultiPoint is given,
    or geometries of type LineString and MultiLineString 
    in the case where 'LineString' or 'MultiLineString' is given.

    Parameters
    --------------------
    vector : SpatialVector
        vector to plot
    ax : matplotlib axis
        the axis to plot the image on
    geomtype : str
        geometry type of the layer
    layer : None | int | str
        layer to plot
    **kwargs (optional)
        any parameter which may be passed to matplotlib.pyplot.plot

    Raises
    --------------------
    ValueError:
        if no layer was specified, and the image contains more than one layer
    ValueError:
        if geomtype is not one of 'Point', 'MultiPoint', 'LineString', 'MultiLineString'
    RuntimeError (from C++):
        if the layer contains a geometry NOT of an acceptable type
    """

    if layer is None:
        if len(vector.layers) == 1:
            layer_name = vector.layers[0]
        else:
            ValueError("no layer was specified, and there is more than one layer in the vector. Specify a layer to plot.");
    elif type(layer) == str:
        layer_name = layer
    elif type(layer) == int:
        layer_name = vector.layers[layer]

    if geomtype == "Point" or geomtype == "MultiPoint":
        points = vector.cpp_vector.get_points(layer_name)
        if 'fmt' in kwargs:
            ax.plot(points[0], points[1], **kwargs)
        else:
            ax.plot(points[0], points[1], '.r', **kwargs) #plot as red points
    elif geomtype == "LineString" or geomtype == "MultiLineString":
        lines = vector.cpp_vector.get_linestrings(layer_name)
        if 'fmt' in kwargs:
            for line in lines:
                ax.plot(line[0], line[1], **kwargs)
        else:
            for line in lines:
                ax.plot(line[0], line[1], '-k', **kwargs) #plot as solid black line
    else:
        raise ValueError("geomtype must be of type 'Point', 'MultiPoint', 'LineString', or 'MultiLineString'");
