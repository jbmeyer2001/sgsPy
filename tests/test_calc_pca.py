import pytest
import numpy as np

import sgs

from files import (
    mraster_geotiff_path,
    pca_result_path
)

class TestPCA:
    #input raster
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    pca_result = sgs.SpatialRaster(pca_result_path)

    @pytest.mark.skip()
    def test_result(self):
        pca = sgs.pca(self.rast, num_comp=3)
        test = pca.band(0)
        correct = self.pca_result.band(0)
        assert np.array_equal(test, correct, equal_nan=True)
        
        test = pca.band(1)
        correct = self.pca_result.band(1)
        assert np.array_equal(test, correct, equal_nan=True)

        test = pca.band(2)
        correct = self.pca_result.band(2)
        assert np.array_equal(test, correct, equal_nan=True)

    @pytest.mark.skip()
    def test_inputs(self):
        pca = sgs.pca(self.rast, num_comp=3)
        pca = sgs.pca(self.rast, num_comp=2)
        pca = sgs.pca(self.rast, num_comp=1)

        with pytest.raises(ValueError):
            pca = sgs.pca(self.rast, num_comp=0)

        with pytest.raises(ValueError):
            pca = sgs.pca(self.rast, num_comp=-1)

        with pytest.raises(ValueError):
            pca = sgs.pca(self.rast, num_comp=4)
