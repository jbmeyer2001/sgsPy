import pytest
import numpy as np

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    sraster_geotiff_path,
    strat_quantiles_zq90_r_path,
    strat_quantiles_pz2_r_path,
)

class TestQuantiles:
    #input raster
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    single_band_rast = sgs.SpatialRaster(sraster_geotiff_path)

    #output raster form running through R version
    zq90_output_rast = sgs.SpatialRaster(strat_quantiles_zq90_r_path)
    pz2_output_rast = sgs.SpatialRaster(strat_quantiles_pz2_r_path)
    
    def test_correct_stratifications_against_R_version(self):
        test_rast = sgs.quantiles(self.rast, num_strata={"zq90": 4})
        test = test_rast.band('strat_zq90')
        correct = np.nan_to_num(np.subtract(self.zq90_output_rast.band(0), 1), nan=-1)
        assert np.array_equal(test, correct, equal_nan=True)

        test_rast = sgs.quantiles(self.rast, num_strata={"pzabove2": [0.2, 0.4, 0.8]})
        test = test_rast.band('strat_pzabove2')
        correct = np.nan_to_num(np.subtract(self.pz2_output_rast.band(0), 1), nan=-1)
        assert np.array_equal(test, correct, equal_nan=True)

    def test_mapping_outputs(self):
        #the python version maps variables differently than the R version
        #so rather than compare against the R version, I'm going to ensure
        #that each mapping stratification corresponds to one single stratification
        #in each of the non-mapping stratifications
        zq90_mapping = {}
        pz2_mapping = {}
        zsd_mapping = {}
        
        test_rast = sgs.quantiles(self.rast, num_strata=[5, [0.2, 0.4, 0.8], 3], map=True)

        for i in range(test_rast.height):
            for j in range(test_rast.width):
                zq90_strat = test_rast.band('strat_zq90')[i, j]
                pz2_strat = test_rast.band('strat_pzabove2')[i, j]
                zsd_strat = test_rast.band('strat_zsd')[i, j]
                map_strat = test_rast.band('strat_map')[i, j]

                if map_strat in zq90_mapping:
                    assert zq90_mapping[map_strat] == zq90_strat
                    assert pz2_mapping[map_strat] == pz2_strat
                    assert zsd_mapping[map_strat] == zsd_strat
                else:
                    zq90_mapping[map_strat] = zq90_strat
                    pz2_mapping[map_strat] = pz2_strat
                    zsd_mapping[map_strat] = zsd_strat

    def test_quantiles_inputs(self):
        test_rast = sgs.quantiles(self.single_band_rast, num_strata=10)
        test_rast = sgs.quantiles(self.single_band_rast, num_strata=[0.00001])

        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(self.single_band_rast, num_strata=[-0.000001, 0.2, 0.4, 0.7])

        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(self.single_band_rast, num_strata=[0.2, 0.4, 0.8, 1.1])

        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(self.rast, num_strata=5)

        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(self.single_band_rast, num_strata=[[0.2, 0.4, 0.8], 5])
        
        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(self.single_band_rast, num_strata=[])
    
    def test_write_functionality(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "rast.tif"
        sgs.quantiles(self.rast, num_strata={'zq90': 4}, filename=str(temp_file))
        test_rast = sgs.SpatialRaster(str(temp_file))
        test = test_rast.band('strat_zq90')
        correct = np.nan_to_num(np.subtract(self.zq90_output_rast.band(0), 1), nan=-1)
        assert np.array_equal(test, correct, equal_nan=True)
