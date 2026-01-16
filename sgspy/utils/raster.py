# ******************************************************************************
#
#  Project: sgs
#  Purpose: GDALDataset wrapper for raster operations
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import importlib.util
import sys
import os
import shutil
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib #for type checking matplotlib.axes.Axes

from .import plot
from .plot import plot_raster

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
print(os.path.join(os.path.dirname(__file__), ".."))
from _sgs import GDALRasterWrapper

#rasterio optional import
try: 
    import rasterio
    RASTERIO = True
except ImportError as e:
    RASTERIO = False

#gdal optional import
try:
    from osgeo import gdal
    from osgeo import gdal_array
    GDAL = True
except ImportError as e:
    GDAL = False

PROJDB_PATH = os.path.join(sys.prefix, "sgspy")

##
# @ingroup user_utils
# This class represents a spatial raster, and is used as an input to many sgs functions. 
# 
# It has a number of additional uses, including accessing the raster data within as a numpy array,
# plotting with matplotlib, as well as converting to a GDAL or Rasterio dataset object. This class
# also has various attributes representing metadata of the raster which may be useful and can be
# seen in the 'Public Attributes' section.
# 
# Accessing raster data
# --------------------
# 
# raster data can be accessed in the form of a NumPy array per band. This can be done using 
# the 'band' function. The band function takes  a single parameter, which must be either
# an integer or a string. If it is an integer, it must refer to a valid zero-indexed band
# number. If it is a string, it must refer to a valid band name within the raster. This function
# may fail if the band is too large to fit in memory.
#
# rast = sgspy.SpatialRaster('test.tif') #raster with three layers
# 
# b0 = rast.band(band=0) @n
# b1 = rast.band(band=1) @n
# b2 = rast.band(band=2) @n
# 
# zq90 = rast.band(band='zq90') @n
# pzabove2 = rast.band(band='pzabove2') @n
# zstd = rast.band(band='zstd') @n
# 
# Accessing raster information
# --------------------
# 
# raster metadata can be displayed using the info() function. Info
# inclues: raster driver, band names, dimensions, pixel size, and bounds.
# 
# rast = sgspy.SpatialRaster('test.tif') @n
# rast.info() 
# 
# Plotting raster
# --------------------
# 
# the plot() function provides a wrapper around matplotlibs imshow 
# functionality (matplotlib.pyplot.imshow). Only a single band can
# be plotted, and for multi-band rasters an indication must be given
# for which band to plot. 
# 
# Target width and heights can be given in the parameters 
# target_width and target_height. Default parameters are 1000 pixels for both. 
# Information on the actual downsampling can be found here:
# https://gdal.org/en/stable/api/gdaldataset_cpp.html#classGDALDataset_1ae66e21b09000133a0f4d99baabf7a0ec
# 
# If no 'band' argument is given, the function will throw an error if the
# image does not contain a single band.
# 
# The 'band' argument allows the end-user to specify either the band
# index or the band name. 'band' may be an int or str.
# 
# Optionally, any of the arguments which may be passed to the matplotlib
# imshow function may also be passed to plot_image(), such as cmap
# for a specific color mapping.
#
# #plots the single band @n
# rast = sgspy.SpatialRaster('test_single_band_raster.tif')  @n
# rast.plot_image()
# 
# #plots the second band @n
# rast = sgspy.SpatialRaster('test_multi_band_raster.tif') @n
# rast.plot(band=1)
# 
# #plots the 'zq90' band @n
# rast = sgspy.SpatialRaster('test_multi_band_raster.tif') @n
# rast.plot(band='zq90')
# 
# Public Attributes
# --------------------
# driver : str @n
#     gdal dataset driver, for info/display purposes @n @n
# width : int @n
#     the pixel width of the raster image @n @n
# height : int @n
#     the pixel height of the raster image @n @n
# band_count : int @n
#     the number of bands in the raster image @n @n
# bands : list[str] @n
#     the raster band names @n @n
# crs : str @n
#     coordinate reference system @n @n
# projection : str @n
#     full projection string as wkt @n @n
# xmin : double @n
#     minimum x value as defined by the gdal geotransform @n @n
# xmax : double @n
#     maximum x value as defined by the gdal geotransform @n @n
# ymin : double @n
#     minimum y value as defined by the gdal geotransform @n @n
# ymax : double @n
#     maximum y value as defined by the gdal geotransform @n @n
# pixel_height : double  @n
#     pixel height as defined by the gdal geotransform @n @n
# pixel_width : double @n
#     pixel width as defined by the gdal geotransform @n @n
# 
# Public Methods
# --------------------
# info() @n
#     takes no arguments, prints raster information to the console @n @n
# plot() @n
#     takes one optional 'band' argument of type int, or str @n @n
# band() @n
#     returns the band data as a numpy array, may throw an error if the raster band is too large
#     
# Optionally, any of the arguments that can be passed to matplotlib.pyplot.imshow 
#     can also be passed to plot_image().
class SpatialRaster:
    
    have_temp_dir = False
    temp_dataset = False
    filename = ""
    closed = False

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
        self.projection
        self.xmin
        self.xmax
        self.ymin
        self.ymax
        self.pixel_height
        self.pixel_width
        self.bands

        Parameters
        --------------------
        image : str
            specifies a raster file path
        """
        if (type(image) is str):
            self.cpp_raster = GDALRasterWrapper(image, PROJDB_PATH)
            self.filename = image
        elif type(image) is GDALRasterWrapper:
            self.cpp_raster = image
        else:
            raise TypeError("'image' parameter of SpatialRaster constructor must be of type str or GDALRasterWrapper")

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
        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

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
        band : str or int
            string representing a band or int representing a band
        """
        if type(band) not in [str, int]:
            raise TypeError("'band' parameter must be of type str or int.")
        
        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

        if type(band) == str:
            band = self.band_name_dict[band]

        return band
 
    def load_arr(self, band_index: int):
        """
        Loads the rasters gdal dataset into a numpy array.
        
        Parameters:
        band : int
            integer representing band index
        """
        if type(band_index) is not int:
            raise TypeError("band_index' parameter must be of type int.")

        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

        self.band_data_dict[band_index] = np.asarray(
            self.cpp_raster.get_raster_as_memoryview(self.width, self.height, band_index).toreadonly(), 
            copy=False
        )

    def band(self, band: str | int):
        """
        gets a numpy array with the specified bands data.

        Parameters:
        band : int | str
            string or int representing band
        """
        if type(band) not in [int, str]:
            raise TypeError("'band' parameter must be of type int or str.")

        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

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
            maximum height in pixels for the image (after downsampling)
        band : int or str
            specification of which bands to plot
        **kwargs
            any parameters which may be passed to matplotlib.pyplot.imshow
        """
        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

        if ax is not None:
            plot_raster(self, ax, target_width, target_width, band, **kwargs)
        else:
            fig, ax = plt.subplots()
            plot_raster(self, ax, target_width, target_width, band, **kwargs)
            plt.show()
       
    @classmethod
    def from_rasterio(cls, ds, arr = None):
        """
        This function is used to convert from a rasterio dataset object representing a raster into an sgspy.SpatialRaster
        object. A np.ndarray may be passed as the 'arr' parameter, if so, the following must be true:
        arr.shape == (ds.count, ds.height, ds.width)

        Examples:

        ds = rasterio.open("rast.tif")
        rast = sgspy.SpatialRaster.from_rasterio(ds)


        ds = rasterio.open("rast.tif")
        arr = ds.read()
        arr[arr < 2] = np.nan
        rast = sgspy.SpatialRaster.from_rasterio(ds, arr)
        """
        if not RASTERIO:
            raise RuntimeError("from_rasterio() can only be called if rasterio was successfully imported, but it wasn't.")

        if type(ds) is not rasterio.io.DatasetReader and type(ds) is not rasterio.io.DatasetWriter:
            raise TypeError("the ds parameter passed to from_raster() must be of type rasterio.io.DatasetReader or rasterio.io.DatasetWriter.")

        if ds.driver == "MEM" and arr is None:
            arr = ds.read()

        if arr is not None:
            if type(arr) is not np.ndarray:
                raise TypeError("the 'arr' parameter, if passed, must be of type np.ndarray")
    
            shape = arr.shape
            if (len(shape)) == 2:
                (height, width) = shape
                if ds.count != 1:
                    raise RuntimeError("if the array parameter contains only a single band with shape (height, width), the raster must contain only a single band.")
            else:
                (band_count, height, width) = shape
                if (band_count != ds.count):
                    raise RuntimeError("the array parameter must contains the same number of bands as the raster with shape (band_count, height, width).")

            if height != ds.height:
                raise RuntimeError("the height of the array passed must be equal to the height of the raster dataset.")

            if width != ds.width:
                raise RuntimeError("the width of the array passed must be equal to the width of the raster dataset.")

            nan = ds.profile["nodata"]
            if nan is None:
                nan = np.nan

            use_arr = True
        else:
            use_arr = False

        if not use_arr:
            #close the rasterio dataset, and open a GDALRasterWrapper of the file
            filename = ds.name
            ds.close()
            return cls(filename)
        else:
            #create an in-memory dataset using the numpy array as the data, and the rasterio dataset to provide metadata
            geotransform = ds.get_transform()
            projection = ds.crs.wkt
            arr = np.ascontiguousarray(arr)
            buffer = memoryview(arr)
            return cls(GDALRasterWrapper(buffer, geotransform, projection, [nan] * ds.count, ds.descriptions, PROJDB_PATH))

    def to_rasterio(self, with_arr = False):
        """
        This function is used to convert an sgspy.SpatialRaster into a rasterio dataset. If with_arr is set to True,
        the function will return a numpy.ndarray as a tuple with the rasterio dataset object.

        Examples:

        rast = sgspy.SpatialRaster('rast.tif')
        ds = rast.to_rasterio()

        rast = sgspy.SpatialRaster('mraster.tif')
        ds, arr = sgs.to_rasterio(with_arr=True)
        """
        if type(with_arr) is not bool:
            raise TypeError("'with_arr' parameter must be of type bool.")

        if not RASTERIO:
            raise RuntimeError("from_rasterio() can only be called if rasterio was successfully imported, but it wasn't.")

        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

        if (self.temp_dataset):
            raise RuntimeError("the dataset has been saved as a temporary file which will be deleted when the C++ object containing it is deleted. the dataset must be either in-memory or have a filename.")

        in_mem = self.driver.find("MEM") != -1

        if with_arr or in_mem:
            bands = []
            for i in range(self.band_count):
                bands.append(np.asarray(self.cpp_raster.get_raster_as_memoryview(self.width, self.height, i)))
            
            #ensure numpy array doesn't accidentally get cleaned up by C++ object deletion
            self.cpp_raster.release_band_buffers()

            arr = np.stack(bands, axis=0)

        if in_mem:
            driver = "MEM"
            width = self.width
            height = self.height
            count = self.band_count
            crs = self.projection
            gt = self.cpp_raster.get_geotransform()
            transform = rasterio.transform.Affine(gt[1], gt[2], gt[0], gt[4], gt[5], gt[3]) #of course the rasterio transform has a different layout than a gdal geotransform... >:(

            dtype = self.cpp_raster.get_data_type()

            if dtype == "":
                raise RuntimeError("sgs dataset has bands with different types, which is not supported by rasterio.")

            nan = self.cpp_raster.get_band_nodata_value(0)

            self.cpp_raster.close()
            self.closed = True

            ds = rasterio.MemoryFile().open(driver=driver, width=width, height=height, count=count, crs=crs, transform=transform, dtype=dtype, nodata=nan)
            ds.write(arr)
        
            for i in range(len(self.bands)):
                ds.set_band_description(i + 1, self.bands[i]) 
        else:
            ds = rasterio.open(self.filename)

        if with_arr:
            return ds, arr
        else:
            return ds

    @classmethod
    def from_gdal(cls, ds, arr = None):
        """
        This function is used to convert from a gdal.Dataset object representing a raster into an sgspy.SpatialRaster
        object. A np.ndarray may be passed as the 'arr' parameter, if so, the following must be true:
        arr.shape == (ds.RasterCount, ds.RasterYSize, ds.RasterXSize)

        Examples:

        ds = gdal.Open("rast.tif")
        rast = sgspy.SpatialRaster.from_gdal(ds)


        ds = gdal.Open("rast.tif")
        bands = []
        for i in range(1, ds.RasterCount + 1):
            bands.append(ds.GetRasterBand(1).ReadAsArray())
        arr = np.stack(bands, axis=0)
        arr[arr < 2] = np.nan
        rast = sgspy.SpatialRaster.from_gdal(ds, arr)
        """
        if not GDAL:
            raise RuntimeError("from_gdal() can only be called if gdal was successfully imported, but it wasn't")

        if type(ds) is not gdal.Dataset:
            raise TypeError("the ds parameter passed to from_gdal() must be of type gdal.Dataset")
    
        if ds.GetDriver().ShortName == "MEM" and arr is None:
            bands = []
            for i in range(1, ds.RasterCount + 1):
                bands.append(ds.GetRasterBand(i).ReadAsArray())
            arr = np.stack(bands, axis=0)

        if arr is not None:
            if type(arr) is not np.ndarray:
                raise TypeError("'arr' parameter, if passed, must be of type np.ndarray")
    
            shape = arr.shape
            if (len(shape)) == 2:
                (height, width) = shape
                if ds.RasterCount != 1:
                    raise RuntimeError("if the array parameter contains only a single band with shape (height, width), the raster must contain only a single band.")
            else:
                (band_count, height, width) = shape
                if (band_count != ds.RasterCount):
                    raise RuntimeError("the array parameter must contains the same number of bands as the raster with shape (band_count, height, width).")

            if height != ds.RasterYSize:
                raise RuntimeError("the height of the array passed must be equal to the height of the raster dataset.")

            if width != ds.RasterXSize:
                raise RuntimeError("the width of the array passed must be equal to the width of the raster dataset.")

            nan_vals = []
            band_names = []
            for i in range(1, ds.RasterCount + 1):
                band = ds.GetRasterBand(i)
                nan_vals.append(band.GetNoDataValue())
                band_names.append(band.GetDescription())

            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()
            arr = np.ascontiguousarray(arr)
            buffer = memoryview(arr)
            
            ds.Close()
            return cls(GDALRasterWrapper(buffer, geotransform, projection, nan_vals, band_names, PROJDB_PATH))
        else:
            filename = ds.GetName()
            
            ds.Close()
            return cls(filename)
            
    def to_gdal(self, with_arr = False):
        """
        This function is used to convert an sgspy.SpatialRaster into a GDAL dataset. If with_arr is set to True,
        the function will return a numpy.ndarray as a tuple with the GDAL dataset object.

        Examples:

        rast = sgspy.SpatialRaster('rast.tif')
        ds = rast.to_gdal()

        rast = sgspy.SpatialRaster('mraster.tif')
        ds, arr = sgs.to_gdal(with_arr=True)
        """
        if not GDAL:
            raise RuntimeError("from_gdal() can only be called if gdal was successfully imported, but it wasn't")

        if self.closed:
            raise RuntimeError("the C++ object which this class wraps has been cleaned up and closed.")

        if (self.temp_dataset):
            raise RuntimeError("the dataset has been saved as a temporary file which will be deleted when the C++ object containing it is deleted. the dataset must be either in-memory or have a filename.")

        in_mem = self.driver.find("MEM") != -1 
        
        if with_arr or in_mem:
            bands = []
            for i in range(self.band_count):
                bands.append(np.asarray(self.cpp_raster.get_raster_as_memoryview(self.width, self.height, i)))
            
            #ensure numpy array doesn't accidentally get cleaned up by C++ object deletion
            self.cpp_raster.release_band_buffers()

            arr = np.stack(bands, axis=0)            

        if in_mem:
            geotransform = self.cpp_raster.get_geotransform()
            projection = self.projection
            nan_vals = []
            for i in range(self.band_count):
                nan_vals.append(self.cpp_raster.get_band_nodata_value(i))
            band_names = self.bands

            self.cpp_raster.close()
            self.closed = True

            ds = gdal.GetDriverByName("MEM").Create("", self.width, self.height, 0, gdal.GDT_Unknown)
            ds.SetGeoTransform(geotransform)
            ds.SetProjection(projection)
            # NOTE: copying from an in-memory GDAL dataset then writing to a rasterio MemoryFile() may cause extra data copying.
            #
            # I'm dong it this way for now instead of somehow passing the data pointer directly, for fear of memory leaks/dangling pointers/accidentally deleting memory still in use.
            for i in range(1, arr.shape[0] + 1):
                band_arr = arr[i - 1]
                ds.AddBand(gdal_array.NumericTypeCodeToGDALTypeCode(band_arr.dtype))
                band = ds.GetRasterBand(i)
                band.WriteArray(band_arr)
                band.SetNoDataValue(nan_vals[i - 1])
                band.SetDescription(band_names[i - 1])
        else:
            self.cpp_raster.close()
            self.closed = True

            ds = gdal.Open(self.filename)

        if with_arr:
            return ds, arr
        else:
            return ds
