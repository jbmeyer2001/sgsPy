import pytest
import numpy as np

import sgs

from files import (
    mraster_geotiff_path,
    inventory_polygons_shapefile_path,
    strat_poly_test1_r_path,
    strat_poly_test2_r_path,
)

class TestPoly:
    #inputs
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    vect = sgs.SpatialVector(inventory_polygons_shapefile_path)

    #output rasters
    test1_output_rast = sgs.SpatialRaster(strat_poly_test1_r_path)
    test2_output_rast = sgs.SpatialRaster(strat_poly_test2_r_path)

    def test_correct_stratifications_against_R_version(self):
        test_rast = sgs.poly(
            self.rast, 
            self.vect, 
            attribute='NUTRIENTS', 
            layer_name='inventory_polygons',
            features=['poor', 'rich', 'medium'],
        )
        test = test_rast[:]
        correct = self.test1_output_rast[:].astype(np.float32)
        correct[correct == 4294967295] = np.nan
        correct = np.subtract(correct, 1)
        assert np.array_equal(test, correct, equal_nan=True)
        
        test_rast = sgs.poly(
            self.rast,
            self.vect,
            attribute='NUTRIENTS',
            layer_name='inventory_polygons',
            features=['poor', ['rich', 'medium']],
        )
        test = test_rast[:]
        correct = self.test2_output_rast[:].astype(np.float32)
        correct[correct == 4294967295] = np.nan
        correct = np.subtract(correct, 1)
        assert np.array_equal(test, correct, equal_nan=True)

    def test_write_functionality(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "rast.tif"
        sgs.poly(
            self.rast, 
            self.vect, 
            attribute='NUTRIENTS', 
            layer_name='inventory_polygons',
            features=['poor', 'rich', 'medium'],
            filename=str(temp_file),
        )
        test_rast = sgs.SpatialRaster(str(temp_file))
        test = test_rast[:]
        correct = self.test1_output_rast[:].astype(np.float32)
        correct[correct == 4294967295] = np.nan
        correct = np.subtract(correct, 1)
        assert np.array_equal(test, correct, equal_nan=True)

