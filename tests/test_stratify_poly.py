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
            'NUTRIENTS', 
            'inventory_polygons',
            ['poor', 'rich', 'medium'],
        )
        for i in range(self.test1_output_rast.height):
            for j in range(self.test1_output_rast.width):
                if np.isnan(test_rast[0, i, j]):
                    #the R image is not a float, it stores nan values as 4294967295
                    assert self.test1_output_rast[0, i, j] == 4294967295
                else:
                    pass
                    #minus 1 to R output because R is 1-indexed
                    assert test_rast[0, i, j] == self.test1_output_rast[0, i, j] - 1
        
        test_rast = sgs.poly(
            self.rast,
            self.vect,
            'NUTRIENTS',
            'inventory_polygons',
            ['poor', ['rich', 'medium']],
        )
        for i in range(self.test2_output_rast.height):
            for j in range(self.test2_output_rast.width):
                if np.isnan(test_rast[0, i, j]):
                    #the R image is not a float, it stores nan values as 4294967295
                    assert self.test2_output_rast[0, i, j] == 4294967295
                else:
                    #minus 1 to R output because R is 1-indexed
                    assert test_rast[0, i, j] == self.test2_output_rast[0, i, j] - 1

    def test_write_functionality(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "rast.tif"
        sgs.poly(
            self.rast, 
            self.vect, 
            'NUTRIENTS', 
            'inventory_polygons',
            ['poor', 'rich', 'medium'],
            filename=str(temp_file),
        )
        test_rast = sgs.SpatialRaster(str(temp_file))
        for i in range(self.test1_output_rast.height):
            for j in range(self.test1_output_rast.width):
                if np.isnan(test_rast[0, i, j]):
                    #the R image is not a float, it stores nan values as 4294967295
                    assert self.test1_output_rast[0, i, j] == 4294967295 
                else:
                    #minus 1 to R output because R is 1-indexed
                    assert test_rast[0, i, j] == self.test1_output_rast[0, i, j] - 1


