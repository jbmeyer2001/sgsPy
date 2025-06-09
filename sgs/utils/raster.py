# ******************************************************************************
#
#  Project: sgs
#  Purpose: GDALDataset wrapper for raster operations
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import json

import numpy as np
import matplotlib.pyplot as plt

from raster import GDALRasterWrapper

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

        the plot_image() function provides a wrapper around matplotlibs imshow 
        functionality (matplotlib.pyplot.imshow). As such, either a single band 
        can be plotted, or three bands can be plotted as an RGB image.

        If no 'bands' argument is given, the function may throw an error if the
        image does not contain 1 or 3 bands.

        The 'bands' argument allows the end-user to specify either the band
        index or the band name for either a scalar or rgb image. 'bands' may
        be an int, str, list, or dict.

        Optionally, any of the arguments which may be passed to the matplotlib
        imshow function may also be passed to plot_image(), such as cmap
        for a specific color mapping.

        examples:
        rast = sgs.utils.raster.SpatialRaster('test.tif')

        #plots either scalar or RGB image depending on number of bands 
        rast.plot_image()

        #plots the second band (index 1) as a scalar image
        rast.plot_image(bands=1)

        #plots the zq90 band as a scalar image
        rast.plot_image(bands='zq90')

        #plots an RGB image with 0 as red, 2 as green, and 1 as blue
        rast.plot_image(bands=[0,2,1])

        #plots an RGB image with red green and blue bands specified
        rast.plot_image(bands={'red':0, 'green':'zsd', 'blue':1}

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
        takes one optional 'bands' argument of type int, str, list, or dict specifying
        the bands to be used in the scalar or RGB image. 
        Optionally, any of the arguments that can be passed to matplotlib.pyplot.imshow 
        can also be passed to plot_image().
    """

    def __init__(self, image):
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
        
        arr is initally set to None, as the array is loaded into a NumPy 
        array only if it is required.

        Parameters
        --------------------
        image : str
            specifies a raster file path

        Raises
        --------------------
        TypeError: 
            if 'image' parameter is not of type str
        RuntimeError (from C++):
            if dataset is not initialized correctly
        RuntimeError (from C++):
            if unable to getgeotransform
        RuntimeError (from C++):
            if unable to get coordinate reference system
        """
        if type(image) == str:
            self.cpp_raster = GDALRasterWrapper(image)
        else:
            raise TypeError(f"SpatialRaster does not accept input of type {type(image)}")

        self.driver = self.cpp_raster.get_driver()
        self.width = self.cpp_raster.get_width()
        self.height = self.cpp_raster.get_height()
        self.band_count = self.cpp_raster.get_band_count()
        self.crs = json.loads(self.cpp_raster.get_crs())
        self.xmin = self.cpp_raster.get_xmin()
        self.xmax = self.cpp_raster.get_xmax()
        self.ymin = self.cpp_raster.get_ymin()
        self.ymax = self.cpp_raster.get_ymax()
        self.pixel_width = self.cpp_raster.get_pixel_width()
        self.pixel_height = self.cpp_raster.get_pixel_height()
        self.arr = None
        self.band_name_dict = {}
        self.bands = self.cpp_raster.get_bands()
        for i in range(0, len(self.bands)):
            self.band_name_dict[self.bands[i]] = i

    def info(self):
        """
        Displays driver, band, size, pixel size, and bound information of the raster.
        """
        print("driver: {}".format(self.driver))
        print("bands: {}".format(*self.bands))
        print("size: {} x {} x {}".format(self.band_count, self.width, self.height))
        print("pixel size: (x, y): ({}, {})".format(self.pixel_height, self.pixel_width))
        print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(self.xmin, self.xmax, self.ymin, self.ymax))

    def load_arr(self):
        """
        Loads the rasters gdal dataset into a numpy array.

        Raises
        --------------------
        RuntimeError (from C++):
            if unable to read raster band
        """
        self.arr = np.asarray(
            self.cpp_raster.get_raster_as_memoryview(self.width, self.height).toreadonly(), 
            copy=False
        )

    def get_band_index(self, band):
        """
        Utilizes the band_name_dict to convert a band name to an index if requried.

        Parameters:
        band: str or int
            string representing a band or int representing a band
        """
        if type(band) == str:
            band = self.band_name_dict[band]

        return band

    def __getitem__(self, index):
        """
        Implements numpy array accesses on self.arr attribute, allowing bands
        to be specified by their name as opposed to index if desired.

        Parameters
        --------------------
        index : int, str, tuple, or slice
            the index of the image desired, allowing bands to be specified as strings

        Raises
        ---------------------
        TypeError: 
            if the index is not of type int, str, tuple, or slice.
        RuntimeError (from C++):
            if unable to read raster band
        """
        if self.arr is None:
            self.load_arr()

        if type(index) == tuple:
            index = (self.get_band_index(index[0]),) + index[1:]
        elif type(index) == slice:
            index = slice(self.get_band_index(index.start), self.get_band_index(index.stop))
        elif type(index) == int or type(index) == str:
            index = self.get_band_index(index)
        else:
            raise TypeError("index must be of type int, str, tuple, or slice")

        return self.arr[index]

    def plot(self, target_width=1000, target_height=1000, bands=None, **kwargs):
        """
        Calls sgs.utils.plot_raster() on self.

        Parameters
        --------------------
        max_width : int
            maximum width in pixels for the image (after downsampling)
        max_height : int
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

        fig, ax = plt.subplots()
        sgs.utils.plot_raster(self, ax, max_width, max_height, bands, **kwargs)
        plt.show()
        
