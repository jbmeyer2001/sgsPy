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

    def test_result(self):
        pca = sgs.pca(self.rast, num_comp=3)
        test = pca.band(0)
        correct = self.pca_result.band(0)
        np.testing.assert_almost_equal(correct, test, decimal=5)
        
        test = pca.band(1)
        correct = self.pca_result.band(1)
        np.testing.assert_almost_equal(correct, test, decimal=5)

        test = pca.band(2)
        correct = self.pca_result.band(2)
        np.testing.assert_almost_equal(correct, test, decimal=5)

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
