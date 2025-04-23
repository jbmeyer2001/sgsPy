import pytest
import sgs
from osgeo import gdal

mraster_small_path = '/home/jbmeyer/extdata/mraster.tif'

class TestSpatialRaster:
    def test_construct_from_path(self):
        rast = sgs.utils.raster.SpatialRaster(mraster_small_path)

    def test_construct_from_gdal_dataset(self):
        dataset = gdal.Open(mraster_small_path, gdal.GA_ReadOnly)
        rast = sgs.utils.raster.SpatialRaster(dataset)

