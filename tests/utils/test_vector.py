import pytest
import sgs

from files import (
    access_shapefile_path,
    existing_shapefile_path,
    existing_geodatabase_path,
    existing_geojson_path,
    inventory_polygons_shapefile_path,
)

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

    @pytest.mark.skip(reason="C++ vector functionality not implemented yet")
    def test_construct_from_path(self):
        vec = sgs.utils.vector.SpatialVector(access_shapefile_path)
        self.access_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(existing_shapefile_path)
        self.existing_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(inventory_polygons_shapefile_path)
        self.inventory_polygons_parameters_check(vec)

    @pytest.mark.skip(reason="C++ vector functionality not implemented yet")
    def test_common_file_formats(self):
        vec = sgs.utils.vector.SpatialVector(existing_geodatabase_path)
        self.existing_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(existing_geojson_path)
        self.existing_parameters_check(vec)
