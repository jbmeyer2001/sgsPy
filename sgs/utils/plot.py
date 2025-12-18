# ******************************************************************************
#
#  Project: sgs
#  Purpose: Plotting rasters and vectors with matplotlib.pyplot
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib #for typing matplotlib.axes.Axes

def plot_raster(raster, 
                ax: matplotlib.axes.Axes, 
                target_width: int = 1000, 
                target_height: int = 1000, 
                band: Optional[int | str] = None, 
                **kwargs):
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
    band (optional) : int or str
        specification of which band to plot
    **kwargs (optional)
        any parameters which may be passed to matplotlib.pyplot.imshow
    """
    #get bands argument as list of int
    if band is None:
        if raster.band_count > 1:
            raise ValueError("'band' argument must be given if raster contains more than one band.")
        band = 0
    else:
        band = raster.get_band_index(band)
    title = raster.bands[band]
    
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
        downsampled_width = int(raster.width / 8)
        downsampled_height = int(raster.height / 8)

    #get the raster data from the cpp object as a numpy array, and ensure no data is nan
    no_data_val = raster.cpp_raster.get_band_nodata_value(band)
    arr = np.asarray(
        raster.cpp_raster.get_raster_as_memoryview(downsampled_width, downsampled_height, band),
        copy=False
    ).astype(np.float64, copy=True)
    arr[arr == no_data_val] = np.nan

    #get raster origin and raster extent
    extent = (raster.xmin, raster.xmax, raster.ymin, raster.ymax) #(left, right, top, bottom)

    #add image to matplotlib
    plt.title(label=title)
    ax.imshow(arr, origin='upper', extent=extent, **kwargs)

def plot_vector(vector, 
                ax: matplotlib.axes.Axes, 
                geomtype: str, 
                layer: Optional[int | str] = None, 
                **kwargs):
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

    if type(layer) == str:
        layer_name = layer
    elif type(layer) == int:
        layer_name = vector.layers[layer]
    elif len(vector.layers) == 1: #layer is None
        layer_name = vector.layers[0]
    else:
        ValueError("no layer was specified, and there is more than one layer in the vector. Specify a layer to plot.");
    
    if geomtype == "Point" or geomtype == "MultiPoint":
        points = vector.cpp_vector.get_points(layer_name)
        if 'fmt' in kwargs:
             ax.plot(points[0], points[1], **kwargs)
        else:
            ax.plot(points[0], points[1], '.r', **kwargs) #specify format as red points if format was not given
    elif geomtype == "LineString" or geomtype == "MultiLineString":
        lines = vector.cpp_vector.get_linestrings(layer_name)
        if 'fmt' in kwargs:
            for line in lines:
                ax.plot(line[0], line[1], **kwargs)
        else:
            for line in lines:
                ax.plot(line[0], line[1], '-k', **kwargs) #specify format as black lines if format was not give
    else:
        raise ValueError("geomtype must be of type 'Point', 'MultiPoint', 'LineString', or 'MultiLineString'");
