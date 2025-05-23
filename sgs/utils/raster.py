import copy
import json
import os

from osgeo import gdal, osr
import numpy as np
import matplotlib.pyplot as plt

class SpatialRaster:
    def __init__(self, image):
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
        for i in range(0, self.layers): #raster bands are 1 indexed
            self.band_name_dict[self.dataset.GetRasterBand(i + 1).GetDescription()] = i

    def info(self):
        print("driver: {}".format(self.driver))
        print("size: {} x {} x {}".format(self.width, self.height, self.layers))
        if self.geotransform:
            print("pixel size: (x, y): ({}, {})".format(self.pixel_height, self.pixel_width))
            print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(self.xmin, self.xmax, self.ymin, self.ymax))

    # NOTE -- this is the most naive implementation possible. There is no manipulation
    # of tiles, and no specification of how much memory to allocate, etc.
    #
    # in the future, for more intensive tasks, this may have to be changed.
    def load_arr(self):
        self.arr = self.dataset.GetVirtualMemArray()

        #array dimensions should always be 3, even for a single band raster
        if self.arr.ndim == 2:
            self.arr = np.array([self.arr])

        if self.arr.ndim != 3:
            raise ValueError("raster image array must have three dimensions.")

    def get_band_index(self, band):
        if type(band) == str:
            band = self.band_name_dict[band]
        elif type(band) == int:
            band = band - 1 #python has arrays 0-indexed, but conventionally bands are 1-indexed

        return band

    def __getitem__(self, index):
        if self.arr is None:
            self.load_arr()

        if type(index) == tuple:
            index = (self.get_band_index(index[0]),) + index[1:]
        elif type(index) == slice:
            index = slice(self.get_band_index(index.start), self.get_band_index(index.stop))
        elif type(index) == int or type(index) == str:
            index = self.get_band_index(index)
        else:
            raise TypeError("__getitem__ index must be of type int, str, tuple, or slice")

        return self.arr[index]

    def arrange_bands_from_list(self, bands_list):
        num_items = len(bands_list)
        if num_items != 1 and num_items != 3:
            raise ValueError("number of bands must be either 1 (for scalar iamges) or 3 (for RGB images).")

        for i in range(num_items):
            bands_list[i] = self.get_band_index(bands_list[i])

        return bands_list

    def arrange_bands_from_dict(self, bands_dict):
        if len(bands_dict) != 3 or 'red' not in bands_dict or 'green' not in bands_dict or 'blue' not in bands_dict:
            raise ValueError("if bands is a dict, it must to have three items with the keys 'red', 'green', and 'blue'.")

        return [
            self.get_band_index(bands_dict["red"]),
            self.get_band_index(bands_dict["green"]),
            self.get_band_index(bands_dict["blue"])
        ]

    def plot_image(self, bands=None, **kwargs):
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
