import pytest
import sgs
from osgeo import gdal

mraster_path = '/home/jbmeyer/extdata/mraster.tif'
mraster_small_path = '/home/jbmeyer/extdata/mraster_small.tif'
sraster_path = '/home/jbmeyer/extdata/sraster.tif'

class TestSpatialRaster:
    def mraster_parameters_check(self, rast):
        assert rast.width == 373
        assert rast.height == 277
        assert rast.layers == 3
        assert rast.pixel_width == 20
        assert rast.pixel_height == 20
        assert rast.xmin == 431100
        assert rast.xmax == 438560
        assert rast.ymin == 5337700
        assert rast.ymax == 5343240

    def mraster_small_parameters_check(self, rast):
        assert rast.width == 141
        assert rast.height == 110
        assert rast.layers == 3
        assert rast.pixel_width == 20
        assert rast.pixel_height == 20
        assert rast.xmin == 432440
        assert rast.xmax == 435260
        assert rast.ymin == 5340000
        assert rast.ymax == 5342200
    
    def sraster_parameters_check(self, rast):
        assert rast.width == 373
        assert rast.height == 277
        assert rast.layers == 1
        assert rast.pixel_width == 20
        assert rast.pixel_height == 20
        assert rast.xmin == 431100
        assert rast.xmax == 438560
        assert rast.ymin == 5337700
        assert rast.ymax == 5343240

    def test_construct_from_path(self):
        rast = sgs.utils.raster.SpatialRaster(mraster_path)
        self.mraster_parameters_check(rast)

        rast = sgs.utils.raster.SpatialRaster(mraster_small_path)
        self.mraster_small_parameters_check(rast)

        rast = sgs.utils.raster.SpatialRaster(sraster_path)
        self.sraster_parameters_check(rast)

    def test_construct_from_gdal_dataset(self):
        dataset = gdal.Open(mraster_path, gdal.GA_ReadOnly)
        rast = sgs.utils.raster.SpatialRaster(dataset)
        self.mraster_parameters_check(rast)

        dataset = gdal.Open(mraster_small_path, gdal.GA_ReadOnly)
        rast = sgs.utils.raster.SpatialRaster(dataset)
        self.mraster_small_parameters_check(rast)

        dataset = gdal.Open(sraster_path, gdal.GA_ReadOnly)
        rast = sgs.utils.raster.SpatialRaster(dataset)
        self.sraster_parameters_check(rast)

'''
TODO: 
 - use test image that has different pixel_width and pixel_height values
 - test the actual pixel values within the image
 - test different image input types
 - test image plotting?
'''
