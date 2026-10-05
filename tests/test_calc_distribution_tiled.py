import numpy as np
import pytest

import sgspy as sgs

from files import mraster_small_tiled_neg_path


class TestDistributionTiledRaster:
    def test_all_negative_band_matches_numpy(self):
        rast = sgs.SpatialRaster(mraster_small_tiled_neg_path)
        bins, counts = sgs.distribution(rast, band=0, bins=5, plot=False)['population']

        values = rast.band(0)
        values = values[~np.isnan(values)]
        assert values.max() <= -1
        assert bins[0] == pytest.approx(values.min())
        assert bins[-1] == pytest.approx(values.max())
        assert sum(counts) == values.size
