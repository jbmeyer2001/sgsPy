# ******************************************************************************
#
#  Project: sgs
#  Purpose: GDALDataset wrapper for vector operations
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

from typing import Optional

import matplotlib.pyplot as plt
import matplotlib #fpr type checking matplotlib.axes.Axes

from.import plot
from .plot import plot_vector

from _sgs import GDALVectorWrapper

class SpatialVector:
    """
    A wrapper class of a GDAL vector dataset.
    This class is primarily used under the hood, although it has a function
    for displaying vector information.

    Accessing vector info:
        
            vector metadata can be displayed using the info() function. All layers
            are displayed unless a specific layer is specified. The per-layer info
            includes: name, number of features, number of fields, geomtype, and bounds.

    Public Attributes:
    --------------------
    layer_names : list[str]
        a list of layer names
    
    Public Methods:
    --------------------
    info()
        takes an optional argument specify the band, and prints vector metadata to console
    """
    def __init__(self, 
                 image: str | GDALVectorWrapper):
        """
        Constructing method for the SpatialVector class.

        Has one required parameter to specify a gdal dataset. The following
        attributes are populated:
        self.cpp_vector
        self.layer_names

        Parameters
        --------------------
        image: str | GDALVectorWrapper
           specifies a path to a vector file or the C++ class object itself

        Raises
        --------------------
        RuntimeError (from C++):
            if dataset is not initialized correctly 
        """
        if (type(image) is str):
            self.cpp_vector = GDALVectorWrapper(image)
        else:
            self.cpp_vector = image

        self.layers = self.cpp_vector.get_layer_names()

    def print_info(self, 
                   layer_name: str, 
                   layer_info: dict):
        """
        prints layer information using the layer_info from self.cpp_vector.

        Parameters
        --------------------
        name : str
            str containing the layer name
        layer_info : dict
            dict containing 'feature_count', 'field_count', 'geometry_type', 'xmax', 'xmin', 'ymax', and 'ymin' items
        """
        print("{} layer info:".format(layer_name))
        print("feature count: {}".format(layer_info['feature_count']))
        print("field count: {}".format(layer_info['field_count']))
        print("geometry type: {}".format(layer_info['geometry_type']))
        print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(
            layer_info['xmin'], 
            layer_info['xmax'], 
            layer_info['ymin'],
            layer_info['ymax']
        ))
        if layer_info['crs']: print("crs: {}".format(layer_info['crs']))
        print()

    def info(self, 
             layer: Optional[int | str] = None):
        """
        calls self.print_info depending on layer parameter. If no layer is given,
        print all layers. A layer may be specified by either a str or an int.

        Parameters
        --------------------
        layer (optional) : str or int
            specifies the layer to print information on
        """
        if type(layer) == str:
            self.print_info(layer, self.cpp_vector.get_layer_info(layer))
        elif type(layer) == int:
            self.print_info(self.layers[layer], self.cpp_vector.get_layer_info(self.layers[layer]))
        else:
            for layer in self.layers:
                self.print_info(layer, self.cpp_vector.get_layer_info(layer))

    def samples_as_wkt(self):
        """
        Calls get_wkt_points on the underlying cpp class, to return
        the samples as wkt strings. 

        This function requires that there be a layer named 'samples' which
        is comprised entirely of Points or MultiPoints. These conditions
        will be satisfied if this SpatialVector is the output of one of the
        sampling functions in the sgs package.

        Raises
        --------------------
        ValueError:
            if this vector does not have a layer called 'samples'
        RuntimeError (from C++):
            if the 'samples' layer has at least one geometry other than Point or MultiPoint
        """
        if "samples" not in self.layers:
            print("this vector does not have a layer 'samples'")
        else:
            return self.cpp_vector.get_wkt_points('samples')

    def plot(self,
        geomtype: str,
        ax: Optional[matplotlib.axes.Axes] = None,
        layer: Optional[int | str] = None, 
        **kwargs):
        """
        Calls plot_vector on self.

        Paramters
        --------------------
        ax : matplotlib.axes.Axes
            axes to plot the raster on
        geomtype : str
            the geometry type to try to print
        layer : None | int | str
            specification of which layer to print
        **kwargs (optional)
            any parameter which may be passed ot matplotlib.pyplot.plot

        Raises
        --------------------
        ValueError:
            if no layer was specified, and the image contains more than one layer
        ValueError:
            if geomtype is not one of 'Point', 'MultiPoint', 'LineString', 'MultiLineString'
        RuntimeError (from C++):
            if the layer contains a geometry NOT of an acceptable type
        """

        if ax is not None: 
            plot_vector(self, ax, geomtype, layer, **kwargs)
        else:
            fig, ax = plt.subplots()
            plot_vector(self, ax, geomtype, layer, **kwargs)
            plt.show()

