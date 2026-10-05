import pytest

import sgspy as sgs

from files import (
    mraster_small_geotiff_path,
    mraster_small_shifted_path,
    nodata_only_path,
)


class TestMapInputValidation:
    def test_rejects_different_geotransform(self):
        a = sgs.breaks(sgs.SpatialRaster(mraster_small_geotiff_path), breaks={'zq90': [3, 5, 11, 18]})
        b = sgs.breaks(sgs.SpatialRaster(mraster_small_shifted_path), breaks={'zsd': [2, 5]})
        with pytest.raises(RuntimeError, match="geotransform"):
            sgs.map((a, 'strat_zq90', 5), (b, 'strat_zsd', 3))

    def test_rejects_zero_threads(self):
        strat = sgs.breaks(sgs.SpatialRaster(mraster_small_geotiff_path), breaks={'zq90': [3, 5, 11, 18]})
        with pytest.raises(ValueError):
            sgs.map((strat, 'strat_zq90', 5), thread_count=0)


class TestQuantilesInputValidation:
    def test_all_nodata_band_raises(self):
        with pytest.raises(RuntimeError, match="only nodata"):
            sgs.quantiles(sgs.SpatialRaster(nodata_only_path), quantiles=4)
