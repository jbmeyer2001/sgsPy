# ******************************************************************************
#
#  Project: sgs
#  Purpose: GDALDataset wrapper for raster operations
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import os
import sys
import platform

if (platform.system() == 'Windows'):
    bin_path = os.path.join(sys.prefix, "sgs")
    os.add_dll_directory(bin_path)

    if bin_path not in os.environ['PATH']:
        os.environ['PROJ_LIB'] = bin_path
 

import json
import shutil
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib #for type checking matplotlib.axes.Axes

from raster import GDALRasterWrapper

from .import plot
from .plot import plot_raster

class SpatialRaster:
    """
    A Python wrapper of the C++ class GDALRasterWrapper. GDAL is used on the C++ side rather
    than the Python side, as it means gdal does not have to be a python dependency. It also
    allows smoother integration with other C++ code.

    This class is primarily used under the hood, although it contains raster metadata, and 
    functions for displaying info and plotting images which are intended for the end-user.

    Accessing raster data:

        raster data can be accessed in the form of a NumPy array using the Python bracket syntax.
        The first dimension is the band, followed by the y and x dimensions like so [band][y][x]. 
        The band can be identified by either a band name or an index, which, following Python's
        standard is zero-indexed. NumPy's array slicing may also be used.

        examples:
        rast = sgs.utils.raster.SpatialRaster('test.tif')

        #returns a 2d numpy array of the first raster band
        rast[0]

        #returns a 3d numpy array of the first and second raster band
        rast[0:2]

        #returns a 2d numpy array of the 'zq90' band
        rast['zq90']

        #returns a 1d numpy array of the first row of the 'zsd' band
        rast['zsd', 0]

        #returns a 1d numpy array of the first column of the 'zsd' band
        rast['zsd', :, 0]

    Accessing raster information:

        raster metadata can be displayed using the info() function. Info
        inclues: raster driver, band names, dimensions, pixel size, and bounds.

        example:
        rast = sgs.utils.raster.SpatialRaster('test.tif')
        rast.info()

    Plotting raster:

        the plot() function provides a wrapper around matplotlibs imshow 
        functionality (matplotlib.pyplot.imshow). Only a single band can
        be plotted, and for multi-band rasters an indication must be given
        for which band to plot. 

        Target width and heights can be given in the parameters 
        target_width and target_height. Default parameters are 1000 pixels for both. 
        Information on the actual downsampling can be found here:
        https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec

        If no 'band' argument is given, the function will throw an error if the
        image does not contain a single.

        The 'band' argument allows the end-user to specify either the band
        index or the band name. 'band' may be an int or str.

        Optionally, any of the arguments which may be passed to the matplotlib
        imshow function may also be passed to plot_image(), such as cmap
        for a specific color mapping.

        examples:
        #plots the single band 
        rast = sgs.SpatialRaster('test_single_band_raster.tif') 
        rast.plot_image()

        #plots the second band
        rast = sgs.SpatialRaster('test_multi_band_raster.tif')
        rast.plot_image(band=1)

        #plots the 'zq90' band
        rast = sgs.SpatialRaster('test_multi_band_raster.tif')
        rast.plot_image(band='zq90')

    Public Attributes
    --------------------
    driver : str
        gdal dataset driver, for info/display purposes
    width : int
        the pixel width of the raster image
    height : int
        the pixel height of the raster image
    band_count : int
        the number of bands in the raster image
    bands : list[str]
        the raster band names
    crs : dict
        dictionary representing coordinate reference system information
    xmin : double
        minimum x value as defined by the gdal geotransform
    xmax : double
        maximum x value as defined by the gdal geotransform
    ymin : double
        minimum y value as defined by the gdal geotransform
    ymax : double
        maximum y value as defined by the gdal geotransform
    pixel_height : double 
        pixel height as defined by the gdal geotransform
    pixel_width : double
        pixel width as defined by the gdal geotransform

    Public Methods
    --------------------
    info()
        takes no arguments, prints raster information to the console
    plot_image()
        takes one optional 'band' argument of type int, or str
        
    Optionally, any of the arguments that can be passed to matplotlib.pyplot.imshow 
        can also be passed to plot_image().
    """
    have_temp_dir = False

    def __init__(self, 
                 image: str | GDALRasterWrapper):
        """
        Constructing method for the SpatialRaster class.

        Has one required parameter to specify a raster path. The following
        attributes are populated:
        self.cpp_raster
        self.driver
        self.width
        self.height
        self.band_count
        self.crs
        self.xmin
        self.xmax
        self.ymin
        self.ymax
        self.pixel_height
        self.pixel_width
        self.arr
        self.bands
        self.band_name_dict
        self.temp_dir
        
        arr is initally set to None, as the array is loaded into a NumPy 
        array only if it is required.

        Parameters
        --------------------
        image : str
            specifies a raster file path

        Raises
        --------------------
        RuntimeError (from C++):
            if dataset is not initialized correctly
        RuntimeError (from C++):
            if unable to getgeotransform
        RuntimeError (from C++):
            if unable to get coordinate reference system
        """
        if (type(image) is str):
            self.cpp_raster = GDALRasterWrapper(image)
        else:
            self.cpp_raster = image

        self.driver = self.cpp_raster.get_driver()
        self.width = self.cpp_raster.get_width()
        self.height = self.cpp_raster.get_height()
        self.band_count = self.cpp_raster.get_band_count()
        self.crs = self.cpp_raster.get_crs()
        self.projection = self.cpp_raster.get_projection().encode('ascii', 'ignore').decode('unicode_escape')
        self.xmin = self.cpp_raster.get_xmin()
        self.xmax = self.cpp_raster.get_xmax()
        self.ymin = self.cpp_raster.get_ymin()
        self.ymax = self.cpp_raster.get_ymax()
        self.pixel_width = self.cpp_raster.get_pixel_width()
        self.pixel_height = self.cpp_raster.get_pixel_height() 
        self.band_name_dict = {}
        self.band_data_dict = {}
        self.bands = self.cpp_raster.get_bands()
        for i in range(0, len(self.bands)):
            self.band_name_dict[self.bands[i]] = i

    def __del__(self):
        if self.have_temp_dir:
            shutil.rmtree(self.temp_dir)

    def info(self):
        """
        Displays driver, band, size, pixel size, and bound information of the raster.
        """
        print("driver: {}".format(self.driver))
        print("bands: {}".format(*self.bands))
        print("size: {} x {} x {}".format(self.band_count, self.width, self.height))
        print("pixel size: (x, y): ({}, {})".format(self.pixel_height, self.pixel_width))
        print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(self.xmin, self.xmax, self.ymin, self.ymax))
        print("crs: {}".format(self.crs))

    def get_band_index(self, band: str | int):
        """
        Utilizes the band_name_dict to convert a band name to an index if requried.

        Parameters:
        band: str or int
            string representing a band or int representing a band
        """
        if type(band) == str:
            band = self.band_name_dict[band]

        return band
 
    def load_arr(self, band_index: int):
        """
        Loads the rasters gdal dataset into a numpy array.

        Raises
        --------------------
        RuntimeError (from C++):
            if unable to read raster band
        RuntimeError (from C++):
            if the band is larger than a gigabyte sgs will not load it into memory
        """

        self.band_data_dict[band_index] = np.asarray(
            self.cpp_raster.get_raster_as_memoryview(self.width, self.height, band_index).toreadonly(), 
            copy=False
        )

    def band(self, band: str | int):
        """
        gets a numpy array with the specified bands data.

        May call load_arr if this data has not been directly
        accessed by Python before.

        Raises
        --------------------
        RuntimeError (from C++):
            if unable to read raster band
        RuntimeError (from C++):
            if the band is larger than a gigabyte sgs will not load it into memory
        """
        index = self.get_band_index(band)

        if index not in self.band_data_dict:
            self.load_arr(index)
        
        return self.band_data_dict[index]

    def plot(self, 
             ax: Optional[matplotlib.axes.Axes] = None,
             target_width: int = 1000, 
             target_height: int = 1000, 
             band: Optional[int | str] = None, 
             **kwargs):
        """
        Calls plot_raster() on self.

        Parameters
        --------------------
        ax : matplotlib.axes.Axes
            axes to plot the raster on
        target_width : int
            maximum width in pixels for the image (after downsampling)
        target_height : int
            maximum height in pixeils for the image (after downsampling)
        band (optional) : int or str
            specification of which bands to plot
        **kwargs (optional)
            any parameters which may be passed to matplotlib.pyplot.imshow

        Raises
        --------------------
        RuntimeError (from C++)
            if unable to read raster band
        """

        if ax is not None:
            plot_raster(self, ax, target_width, target_width, band, **kwargs)
        else:
            fig, ax = plt.subplots()
            plot_raster(self, ax, target_width, target_width, band, **kwargs)
            plt.show()
        
