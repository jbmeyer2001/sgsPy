# ******************************************************************************
#
#  Project: sgs
#  Purpose: GDALDataset wrapper for vector operations
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import sys
import os
import tempfile
from typing import Optional
import warnings

import matplotlib.pyplot as plt
import matplotlib #fpr type checking matplotlib.axes.Axes

from.import plot
from .plot import plot_vector

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from _sgs import GDALVectorWrapper

try:
    import geopandas as gpd
    GEOPANDAS = True
except ImportError as e:
    GEOPANDAS = False

PROJDB_PATH = os.path.join(sys.prefix, "sgspy")

##
# @ingroup user_utils
# This class represents a spatial vector, and is used as an input to many sgs functions.
#
# It has a number of additional uses, including displaying info about the vector, converting 
# to a GDAL or GeoPandas object. 
# 
# Accessing vector info:
#     
# vector metadata can be displayed using the info() function. All layers
# are displayed unless a specific layer is specified. The per-layer info
# includes: name, number of features, number of fields, geomtype, and bounds.
# 
# Public Attributes:
# --------------------
# layer_names : list[str] @n
#     a list of layer names
# 
# Public Methods:
# --------------------
# info() @n
#     takes an optional argument specify the band, and prints vector metadata to console
class SpatialVector:

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
        """
        if type(image) is str:
            self.cpp_vector = GDALVectorWrapper(image, PROJDB_PATH)
        elif type(image) is GDALVectorWrapper:
            self.cpp_vector = image
        else:
            raise TypeError("'image' parameter to SpatialVector constructor must be of type str or GDALVectorWrapper.")

        self.layers = self.cpp_vector.get_layer_names()

    def print_info(self, 
                   layer_name: str, 
                   layer_info: dict):
        """
        prints layer information using the layer_info from self.cpp_vector.

        This is an internal function not meant to be used by the end user.

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
        layer : str or int
            specifies the layer to print information on
        """
        if layer is not None and type(layer) not in [int, str]:
            raise TypeError("'layer' parameter, if given, must be of type int or str.")

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
        **kwargs
            any parameter which may be passed ot matplotlib.pyplot.plot
        """

        if ax is not None: 
            plot_vector(self, ax, geomtype, layer, **kwargs)
        else:
            fig, ax = plt.subplots()
            plot_vector(self, ax, geomtype, layer, **kwargs)
            plt.show()

    @classmethod
    def from_geopandas(cls, obj, layer_name: str=None):
        """
        This function is used to convert a geopandas object into an sgspy.SpatialVector. The geopandas object
        may either by of type GeoDataFrame or GeoSeries.

        If a particular layer name is desired, it can be passed as a parameter.

        Examples:

        gdf = gpd.read_file("access.shp")
        access = sgspy.SpatialVector.from_geopandas(gdf)

        
        gs = gpd['geometry'] #geometry column is a geoseries
        access = sgspy.SpatialVector.from_geopandas(gs)


        gdf = gpd.read_file("access.shp")
        gdf = gdf[gdf == "LineString"]
        access = sgspy.SpatialVector.from_geopandas(gdf)
        """
        if layer_name is not None and type(layer_name) is not str:
            raise TypeError("layer_name, if given, must be of type 'str'.")

        if not GEOPANDAS:
            raise RuntimeError("from_geopandas() can only be called if geopandas was successfully imported, but it wasn't.")

        if type(obj) is not gpd.geodataframe.GeoDataFrame and type(obj) is not gpd.geoseries.GeoSeries:
            raise TypeError("the object passed must be of type geopandas GeoDataFrame or GeoSeries")

        #get a the geopandas object as a geoseries
        if type(obj) is gpd.geodataframe.GeoDataFrame:
            if 'geometry' not in obj.columns:
                raise RuntimeError("'geometry' must be a column in the geodataframe passed")
            gdf = obj
            gs = gdf['geometry']
        else:
            gs = obj
            gdf = gpd.GeoDataFrame(gs)
 
        if layer_name is None: layer_name = 'geopandas_geodataframe'
        projection = gs.crs.to_wkt()
       
        # the conversion to geojson may raise a warning about how the projection is unable to be converted to the new format correctly
        #
        # however, we have the projection from the geodataframes geoseries crs, so it won't be an issue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            geojson = gdf.to_json().encode('utf-8')

        cpp_vector = GDALVectorWrapper(geojson, projection, layer_name, PROJDB_PATH)
        return cls(cpp_vector)

    def to_geopandas(self):
        """
        This function is used to convert an sgspy.SpatialVector into a geopandas geodataframe.

        Examples:

        access = sgspy.SpatialVector("access.shp")
        gdf = access.to_geopandas()
        """
        if not GEOPANDAS:
            raise RuntimeError("to_geopandas() can  only be called if geopandas was successfully imported, but it wasn't.")

        tempdir = tempfile.gettempdir()
        file = os.path.join(tempdir, "temp.geojson")

        #get the projection info
        projection = self.cpp_vector.get_projection()

        #write the dataset to a tempfile
        self.cpp_vector.write(file)
    
        # This method of writing to a file then reading from that file is definitely clunky,
        # however it's easy. Theres the possiblity of iterating through every field within
        # every feature, and needing to then call a different function depending on the data
        # type of the field (because C++ types are rigid). That may still be done in the future,
        # but for now this works.

        #have geopandas read the tempfile
        gdf = gpd.read_file(file)
        if projection != "": gdf.set_crs(projection, inplace=True, allow_override=True)

        #remove the geojson file
        os.remove(file)

        return gdf
        
