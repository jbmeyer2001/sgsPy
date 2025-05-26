import copy
import json
import os

from osgeo import gdal, osr
import numpy as np
import matplotlib.pyplot as plt

class SpatialRaster:
    """
    A wrapper class of a GDAL raster dataset. 
    This class is primarily used under the hood, although it has functions for plotting
    and displaying raster info which are intended for the end-user.

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
        gdal datast driver, for info/display purposes
    width : int
        the pixel width of the raster image
    height : int
        the pixel height of the raster image
    layers : int
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

        Has one required parameter to specify a gdal dataset. The following
        attributes are populated using the given dataset:
        self.dataset
        self.driver
        self.width
        self.height
        self.layers
        self.crs
        self.geotransform
        self.xmin
        self.xmax
        self.ymin
        self.ymax
        self.pixel_height
        self.pixel_width
        self.arr
        self.bands
        self.bands_name_dict

        xmin, xmax, ymin, ymax, pixel_height, and pixel_width are only 
        set if geotransform exists.
        
        arr is initally set to None, as the array is loaded into memory 
        only if it is required.

        Parameters
        --------------------
        image : str OR gdal.Dataset
            specifies either a gdal dataset, or a path to a gdal dataset

        Raises
        --------------------
        TypeError: 
            if 'image' parameter is not of type str or gdal.Dataset
        ValueError:
            if dataset is not loaded
        """
        if type(image) == str:
            self.dataset = gdal.Open(image, gdal.GA_ReadOnly)
        elif type(image) == gdal.Dataset:
            self.dataset = image
        else:
            raise TypeError(f"SpatialRaster does not except input of type {type(image)}")

        if not self.dataset:
            raise ValueError("dataset must exist")

        self.driver = f"{self.dataset.GetDriver().ShortName}/{self.dataset.GetDriver().LongName}"
        self.width = self.dataset.RasterXSize
        self.height = self.dataset.RasterYSize
        self.layers = self.dataset.RasterCount
        self.crs = json.loads(osr.SpatialReference(wkt=self.dataset.GetProjection()).ExportToPROJJSON())
        self.geotransform = self.dataset.GetGeoTransform()
        if self.geotransform:
            xbounds = [self.geotransform[0], self.geotransform[0] + self.geotransform[1] * self.width + self.geotransform[2] * self.height]
            ybounds = [self.geotransform[3], self.geotransform[3] + self.geotransform[4] * self.width + self.geotransform[5] * self.height]
            self.xmin = np.min(xbounds)
            self.xmax = np.max(xbounds)
            self.ymin = np.min(ybounds)
            self.ymax = np.max(ybounds)
            self.pixel_height = np.abs(self.geotransform[1])
            self.pixel_width = np.abs(self.geotransform[5])
        self.arr = None
        self.band_name_dict = {}
        self.bands = []
        for i in range(0, self.layers): #raster bands are 1 indexed
            band_name = self.dataset.GetRasterBand(i + 1).GetDescription()
            self.band_name_dict[band_name] = i
            self.bands.append(band_name)

    def info(self):
        """
        Displays driver, band, size, pixel size, and bound information of the dataset.
        """
        print("driver: {}".format(self.driver))
        print("bands: {}".format(*self.bands))
        print("size: {} x {} x {}".format(self.layers, self.width, self.height))
        if self.geotransform:
            print("pixel size: (x, y): ({}, {})".format(self.pixel_height, self.pixel_width))
            print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(self.xmin, self.xmax, self.ymin, self.ymax))

    # TODO: add less naive implementation for resource-intensive tasks (tiling?, etc.)
    def load_arr(self):
        """
        Loads the gdal datsaet into a NumPy array using GetVirtualMemArray().

        Raises
        --------------------
        ValueError: 
            if array does not have three dimensions [band][y][x]
        """
        self.arr = self.dataset.GetVirtualMemArray()

        #array dimensions should always be 3, even for a single band raster
        if self.arr.ndim == 2:
            self.arr = np.array([self.arr])

        if self.arr.ndim != 3:
            raise ValueError("raster image array must have three dimensions.")

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
            the index of the image desired, allowign bands to be specified as strings

        Raises
        ---------------------
        TypeError: 
            if the index is not of type int, str, tuple, or slice.
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

    def arrange_bands_from_list(self, bands_list):
        """
        Used by plot_image function. Converts all bands in initial list to int indexes
        using self.get_band_index().

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
            bands_list[i] = self.get_band_index(bands_list[i])

        return bands_list

    def arrange_bands_from_dict(self, bands_dict):
        """
        Used by plot_image function. Converts a dict which specifies RGB values
        into a list of  indexes using self.get_band_index().

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
            self.get_band_index(bands_dict["red"]),
            self.get_band_index(bands_dict["green"]),
            self.get_band_index(bands_dict["blue"])
        ]

    def plot_image(self, bands=None, **kwargs):
        """
        Plots the specified bands using matplotlib.pyplot.imshow function.

        Parameters
        --------------------
        bands (optional) : int or str or list or dict
            specification of which bands to plot
        **kwargs (optional)
            any parameters which may be passed to matplotlib.pyplot.imshow

        Raises
        --------------------
        TypeError:
            if 'bands' is not of type int, str, list, or dict
        """
        if bands is None:
            bands = self.arrange_bands_from_list([*range(self.layers)])
        elif type(bands) == list:
            bands = self.arrange_bands_from_list(bands)
        elif type(bands) == dict:
            bands = self.arrange_bands_from_dict(bands)
        elif type(bands) in [str, int]: 
            bands = [self.get_band_index(bands)]
        else:
            raise TypeError("'bands' parameter must be of type None, list, dict, str, or int.")

        if self.arr is None:
            self.load_arr()

        display_arr = np.moveaxis(self.arr[bands, :, :], 0, 2)
        plt.imshow(display_arr, **kwargs)
        plt.show()
