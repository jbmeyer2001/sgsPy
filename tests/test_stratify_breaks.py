import pytest
import numpy as np

import sgs

from files import (
    mraster_geotiff_path,
    sraster_geotiff_path,
    strat_breaks_zq90_r_path,
    strat_breaks_pz2_r_path,
)

class TestBreaks:
    #input raster
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    single_band_rast = sgs.SpatialRaster(sraster_geotiff_path)

    #output rasters from running through R version
    zq90_output_rast = sgs.SpatialRaster(strat_breaks_zq90_r_path)
    pz2_output_rast = sgs.SpatialRaster(strat_breaks_pz2_r_path)

    def test_correct_stratifications_against_R_version(self):
        test_rast = sgs.breaks(self.rast, breaks={'zq90': [3, 5, 11, 18]})
        test = test_rast[:]
        correct = np.subtract(self.zq90_output_rast[:], 1)
        assert np.array_equal(test, correct, equal_nan=True)

        test_rast = sgs.breaks(self.rast, breaks={'pzabove2': [20, 40, 60, 80]}) #pzabove2
        test = test_rast[:]
        correct = np.subtract(self.pz2_output_rast[:], 1)
        assert np.array_equal(test, correct, equal_nan=True)
            
    def test_mapping_outputs(self):
        #the python version maps variables differently than the R version
        #so rather than compare against the R version, I'm going to ensure
        #that each mapping stratification corresponds to one single stratification
        #in each of the non-mapping stratifications
        zq90_mapping = {}
        pz2_mapping = {}
        zsd_mapping = {}
        
        test_rast = sgs.breaks(self.rast, breaks=[[3, 5, 11, 18], [40, 60, 80], [2, 5]], map=True)
        
        for i in range(test_rast.height):
            for j in range(test_rast.width):
                zq90_strat = test_rast['strat_zq90', i, j]
                pz2_strat = test_rast['strat_pzabove2', i, j]
                zsd_strat = test_rast['strat_zsd', i, j]
                map_strat = test_rast['strat_map', i, j]

                if map_strat in zq90_mapping:
                    assert zq90_mapping[map_strat] == zq90_strat
                    assert pz2_mapping[map_strat] == pz2_strat
                    assert zsd_mapping[map_strat] == zsd_strat
                else:
                    zq90_mapping[map_strat] = zq90_strat
                    pz2_mapping[map_strat] = pz2_strat
                    zsd_mapping[map_strat] = zsd_strat

    def test_breaks_inputs(self):
        #ensuring no errors with single band input
        test_rast = sgs.breaks(self.single_band_rast, breaks=[1,3])
        test_rast = sgs.breaks(self.single_band_rast, breaks=[1])

        with pytest.raises(ValueError):
            test_rast = sgs.breaks(self.rast, breaks=[2, 5])

        with pytest.raises(ValueError):
            test_rast = sgs.breaks(self.single_band_rast, breaks=[[3, 5, 11, 18], [2, 5]])
        
        with pytest.raises(ValueError):
            test_rast = sgs.breaks(self.single_band_rast, breaks=[])

    def test_write_functionality(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "rast.tif"
        sgs.breaks(self.rast, breaks={'zq90': [3, 5, 11, 18]}, filename=str(temp_file))
        test_rast = sgs.SpatialRaster(str(temp_file))
        test = test_rast[:]
        correct = np.subtract(self.zq90_output_rast[:], 1)
        assert np.array_equal(test, correct, equal_nan=True)
