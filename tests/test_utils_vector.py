import pytest
import sgspy as sgs

from files import (
    access_shapefile_path,
    existing_shapefile_path,
    existing_geodatabase_path,
    existing_geojson_path,
    inventory_polygons_shapefile_path,
)

class TestSpatialVector:
    def access_parameters_check(self, vec):
        assert len(vec.layers) == 1
        assert vec.layers[0] == 'access'
        layer_info = vec.cpp_vector.get_layer_info('access')
        assert int(layer_info["feature_count"]) == 167
        assert int(layer_info["field_count"]) == 2
        assert layer_info["geometry_type"] == "Line String"
        assert float(layer_info["xmin"]) == 431100
        assert float(layer_info["xmax"]) == 438560
        assert float(layer_info["ymin"]) == 5337700
        assert float(layer_info["ymax"]) == 5343240

    def existing_parameters_check(self, vec):
        assert len(vec.layers) == 1
        assert vec.layers[0] == 'existing'
        layer_info = vec.cpp_vector.get_layer_info('existing')
        assert int(layer_info["feature_count"]) == 200
        assert int(layer_info["field_count"]) == 1
        assert layer_info["geometry_type"] == "Point"
        assert float(layer_info["xmin"]) == 431110
        assert float(layer_info["xmax"]) == 438530
        assert float(layer_info["ymin"]) == 5337710
        assert float(layer_info["ymax"]) == 5343230

    def inventory_polygons_parameters_check(self, vec):
        assert len(vec.layers) == 1
        assert vec.layers[0] == 'inventory_polygons'
        layer_info = vec.cpp_vector.get_layer_info('inventory_polygons')
        assert int(layer_info["feature_count"]) == 632
        assert int(layer_info["field_count"]) == 3
        assert layer_info["geometry_type"] == "Polygon"
        assert float(layer_info["xmin"]) == 431100
        assert float(layer_info["xmax"]) == 438560
        assert float(layer_info["ymin"]) == 5337700
        assert float(layer_info["ymax"]) == 5343240

    def test_construct_from_path(self):
        vec = sgs.utils.vector.SpatialVector(access_shapefile_path)
        self.access_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(existing_shapefile_path)
        self.existing_parameters_check(vec)

        vec = sgs.utils.vector.SpatialVector(inventory_polygons_shapefile_path)
        self.inventory_polygons_parameters_check(vec)

    def test_geojson_format(self):
        vec = sgs.utils.vector.SpatialVector(existing_geojson_path)
        self.existing_parameters_check(vec)
