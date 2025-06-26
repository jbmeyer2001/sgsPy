import pytest
import numpy as np

import sgs

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
        test_rast = sgs.quantiles(rast, num_strata={"zq90": 4})
        for i in range(zq90_output_rast.height):
            for j in range(zq90_output_rast.width):
                if np.isnan(test_rast[i][j]):
                    assert np.isnan(zq90_output_rast[i][j])
                else:
                    #plus one to R output because R is 1-indexed
                    assert test_rast[i][j] == zq90_output_rast[i][j] + 1

        test_rast = sgs.quantiles(rast, num_strata={"pzabove2": [0.2, 0.4, 0.8]})
        for i in range(pz2_output_rast.height):
            for j in range(pz2_output_rast.width):
                if np.isnan(test_rast[i][j]):
                    assert np.isnan(pz2_output_rast[i][j])
                else:
                    #plus one to R output because R is 1-indexed
                    assert test_rast[i][j] == zq90_output_rast[i][j] + 1

    def test_mapping_outputs(self):
        #the python version maps variables differently than the R version
        #so rather than compare against the R version, I'm going to ensure
        #that each mapping stratification corresponds to one single stratification
        #in each of the non-mapping stratifications
        zq90_mapping = {}
        pz2_mapping = {}
        zsd_mapping = {}
        
        test_rast = sgs.quantiles(rast, num_strata=[5, [0.2, 0.4, 0.8], 3], map=True)

        for i in range(test_rast.height):
            for j in range(test_rast.width):
                zq90_strat = test_rast['strat_zq90', i, j]
                pz2_strat = test_rast['strat_pzabove2', i, j]
                zsd_strat = test_rast['strat_zsd', i, j]
                map_strat = test_rast['strat_map', i, j]

                if map_strat in zq90_mapping:
                    assert zq90_mapping[map_strat] == zq90_strat
                    assert pz2_mapping[map_strat] == pz2_strat
                    assert zsd_mappingt[map_strat] == zsd_strat
                else:
                    zq90_mapping[map_strat] = zq90_strat
                    pz2_mapping[map_strat] = pz2_strat
                    zsd_mapping[map_strat] = zsd_strat

    def test_quantiles_inputs(self):
        test_rast = sgs.quantiles(single_band_rast, num_strata=10)
        test_rast = sgs.quantiles(single_band_rast, num_strata=[0.00001])

        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(single_band_rast, num_strata=[-0.000001, 0.2, 0.4, 0.7])

        with pytest.raises(ValueError):
            test_rast = sgs.quantiles(single_band_rast, num_strata=[0.2, 0.4, 0.8, 1.1])

        with pytest.expect(ValueError):
            test_rast = sgs.quantiles(rast, num_strata=5)

        with pytest.expect(ValueError):
            test_rast = sgs.quantiles(single_band_rast, num_strata=[[0.2, 0.4, 0.8], 5])
        
        with pytest.expect(ValueError):
            test_rast = sgs.quantiles(single_band_rast, num_strata=[])

    def test_write_functionality(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "rast.tif"
        sgs.breaks(rast, breaks={'zq90': [3, 5, 11, 18]}, filename=str(temp_file))
        test_rast = sgs.SpatialRaster(str(temp_file))
        for i in range(test_rast.height):
            for j in range(test_rast.width):
                if np.isnan(test_rast[i][j]):
                    assert np.isnan(zq90_output_rast[i][j])
                else:
                    #plus one to R output because R is 1-indexed
                    assert test_rast[i][j] == zq90_output_rast[i][j] + 1

