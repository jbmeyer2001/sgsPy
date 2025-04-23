import os

from osgeo import gdal, osr
import numpy as np
import matplotlib.pyplot as plt

class SpatialRaster:
    def __init__(self, image):
        # TODO add other input types as required
        if type(image) == str:
            self.dataset = gdal.Open(image, gdal.GA_ReadOnly)
        elif type(image) == gdal.Dataset:
            self.dataset = image
        elif type(image) == SpatialRaster:
            self.dataset = image.dataset
        else:
            raise TypeError(f"SpatialRaster does not except input of type {type(image)}")

        if not self.dataset:
            raise ValueError("dataset must exist")

        self.driver = f"{self.dataset.GetDriver().ShortName}/{self.dataset.GetDriver().LongName}"
        self.width = self.dataset.RasterXSize
        self.height = self.dataset.RasterYSize
        self.layers = self.dataset.RasterCount
        self.crs = osr.SpatialReference(wkt=self.dataset.GetProjection())
        self.geotransform = self.dataset.GetGeoTransform()
        if self.geotransform:
            xbounds = [self.geotransform[0], self.geotransform[0] + self.geotransform[1] * self.width + self.geotransform[2] * self.height]
            ybounds = [self.geotransform[3], self.geotransform[3] + self.geotransform[4] * self.width + self.geotransform[5] * self.height]
            self.xmin = np.min(xbounds)
            self.xmax = np.max(xbounds)
            self.ymin = np.min(ybounds)
            self.ymax = np.max(ybounds)
        self.arr = None
        self.band_name_dict = {}
        for i in range(0, self.layers): #raster bands are 1 indexed
            self.band_name_dict[self.dataset.GetRasterBand(i + 1).GetDescription()] = i


    def info(self):
        print("driver: {}".format(self.driver))
        print("size: {} x {} x {}".format(self.width, self.height, self.layers))
        if self.geotransform:
            print("pixel size: (x, y): ({}, {})".format(np.abs(self.geotransform[1]), np.abs(self.geotransform[5])))
            print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(self.xmin, self.xmax, self.ymin, self.ymax))

    # NOTE -- this is the most naive implementation possible. There is no manipulation
    # of tiles, and no specification of how much memory to allocate, etc.
    #
    # TODO in the future, for more intensive tasks, this may have to be changed.
    def load_arr(self):
        self.arr = self.dataset.GetVirtualMemArray()

    def get_band_index(self, band):
        if type(band) == str:
            band = self.band_name_dict[band]
        return band

    def __getitem__(self, index):
        if self.arr is None:
            self.load_arr()

        if type(index) == tuple:
            index = (self.get_band_index(index[0]),) + index[1:]
        else:
            index = self.get_band_index(index)

        return self.arr[index]

    def arrange_bands_from_list(bands_list):
        #TODO implement this
        #TODO check to ensure alpha isn't included because I won't use it
        pass

    def arrange_bands_from_dict(bands_dict):
        #TODO implement this
        #only accept lists of length 1 or 3
        pass

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
            pass
            #TODO raise value error

        if self.arr is None:
            self.load_arr()

        display_arr = np.moveaxis(self.arr[bands, :, :], 0, 2)
        plt.imshow(display_arr, **kwargs)
        plt.show()