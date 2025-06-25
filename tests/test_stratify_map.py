import pytest

import sgs

from files import mraster_geotiff_path

class TestMap:
    rast = sgs.SpatialRaster(mraster_geotiff_path)
   
    def test_correct_outputs(self):
        zq90_mapping = {}
        pz2_mapping = {}
        zsd_mapping = {}

        breaks = sgs.breaks(rast, breaks={'zq90': [3, 5, 11, 18], 'pzabove2': [20, 40, 60, 80]})
        quantiles = sgs.quantiles(rast, num_strata={'zsd': 25})
        zq90_rast = breaks['strat_zq90']
        pz2_rast = breaks['strat_pzabove2']
        zsd_rast = quantiles['strat_zsd']

        mapped = sgs.map((breaks, ['strat_zq90', 'strat_pzabove2'], [5, 5]), (quantiles, 'strat_zsd', 25))
        
        assert mapped.height == breaks.height
        assert mapped.width == breaks.width
        for i in range(mapped.height):
            for j in range(mapped.width):
                sq90_strat = zq90_rast[i, j]
                pz2_strat = pz2_rast[i, j]
                zsd_strat = zsd_rast[i, j]
                map_strat = ['strat_map', i, j]

                if map_strat in zq90_mapping:
                    assert zq90_mapping[map_strat] == zq90_strat
                    assert pz2_mapping[map_strat] == pz2_strat
                    assert zsd_mappingt[map_strat] == zsd_strat
                else:
                    zq90_mapping[map_strat] = zq90_strat
                    pz2_mapping[map_strat] = pz2_strat
                    zsd_mapping[map_strat] = zsd_strat


    def test_inputs(self):
        breaks = sgs.breaks(rast, breaks={'zq90': [3, 5, 11, 18], 'pzabove2': [20, 40, 60, 80]})
        quantiles = sgs.quantiles(rast, num_strata={'zsd': 25})
        
        #no error when passed as ints
        mapped = sgs.map((breaks, [0, 1], [5, 5]))

        with pytest.raises(TypeError):
            mapped = sgs.map((breaks, ['strat_zq90', 'strat_pzabove2'], 5))
        
        with pytest.raises(TypeError):
            mapped = sgs.map((breaks, 'strat_zq90', [5]))
        
        with pytest.raises(ValueError):
            mapped = sgs.map((breaks, ['strat_zq90', 'strat_pzabove2'], [5]))
        
        with pytest.raises(ValueError):
            mapped = sgs.map((breaks, ['strat_zq90', 'strat_pzabove2', 'strat_zsd'], [5, 5, 25])

        with pytest.raises(ValueError):
            mapped = sgs.map((breaks, ['strat_zq90', 'pzabove2'], [5, 5]))
     
        with pytest.raises(ValueError):
            mapped = sgs.map((breaks, [1, 2], [5, 5]))

    def test_write_functionality(self, tmp_path):
        zq90_mapping = {}
        pz2_mapping = {}
        zsd_mapping = {}

        breaks = sgs.breaks(rast, breaks={'zq90': [3, 5, 11, 18], 'pzabove2': [20, 40, 60, 80]})
        quantiles = sgs.quantiles(rast, num_strata={'zsd': 25})
        zq90_rast = breaks['strat_zq90']
        pz2_rast = breaks['strat_pzabove2']
        zsd_rast = quantiles['strat_zsd']

        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "rast.tif"
        sgs.map((breaks, ['strat_zq90', 'strat_pzabove2'], [5, 5]), (quantiles, 'strat_zsd', 25), filename=str(temp_file))
        mapped = sgs.SpatialRaster(str(temp_file))

        assert mapped.height == breaks.height
        assert mapped.width == breaks.width
        for i in range(mapped.height):
            for j in range(mapped.width):
                sq90_strat = zq90_rast[i, j]
                pz2_strat = pz2_rast[i, j]
                zsd_strat = zsd_rast[i, j]
                map_strat = ['strat_map', i, j]

                if map_strat in zq90_mapping:
                    assert zq90_mapping[map_strat] == zq90_strat
                    assert pz2_mapping[map_strat] == pz2_strat
                    assert zsd_mappingt[map_strat] == zsd_strat
                else:
                    zq90_mapping[map_strat] = zq90_strat
                    pz2_mapping[map_strat] = pz2_strat
                    zsd_mapping[map_strat] = zsd_strat



