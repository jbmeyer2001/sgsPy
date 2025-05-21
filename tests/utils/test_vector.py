import pytest
import sgs
from osgeo import gdal

access_path = '/home/jbmeyer/extdata/access.shp'
existing_path = '/home/jbmeyer/extdata/existing.shp'
inventory_polygons_path = '/home/jbmeyer/extdata/inventory_polygons.shp'

class TestSpatialVector:
    def access_parameters_check(self, vec):
        info = vec.layers_info[0]
        assert info['name'] == 'access'
        assert info['features'] == 167
        assert info['fields'] == 2
        (xmin, xmax, ymin, ymax) = info['extent']
        assert xmin == 431100
        assert xmax == 438560
        assert ymin == 5337700
        assert ymax == 5343240

    def existing_parameters_check(self, vec):
        info = vec.layers_info[0]
        assert info['name'] == 'existing'
        assert info['features'] == 200
        assert info['fields'] == 1
        (xmin, xmax, ymin, ymax) = info['extent']
        assert xmin == 431110
        assert xmax == 438530
        assert ymin == 5337710
        assert ymax == 5343230

    def inventory_polygons_parameters_check(self, vec):
        info = vec.layers_info[0]
        assert info['name'] == 'inventory_polygons'
        assert info['features'] == 632
        assert info['fields'] == 3
        (xmin, xmax, ymin, ymax) = info['extent']
        assert xmin == 431100
        assert xmax == 438560
        assert ymin == 5337700
        assert ymax == 5343240

    def test_construct_from_path(self):
        vec = sgs.utils.vector.SpatialVector(access_path)
        self.access_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(existing_path)
        self.existing_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(inventory_polygons_path)
        self.inventory_polygons_parameters_check(vec)

    def test_construct_from_gdal_dataset(self):
        dataset = gdal.OpenEx(access_path, gdal.OF_VECTOR)
        vec = sgs.utils.vector.SpatialVector(dataset)
        self.access_parameters_check(vec)

        dataset = gdal.OpenEx(existing_path, gdal.OF_VECTOR)
        vec = sgs.utils.vector.SpatialVector(dataset)
        self.existing_parameters_check(vec)

        dataset = gdal.OpenEx(inventory_polygons_path, gdal.OF_VECTOR)
        vec = sgs.utils.vector.SpatialVector(dataset)
        self.inventory_polygons_parameters_check(vec)
