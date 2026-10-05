import pytest

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    sraster_geotiff_path,
)


class TestStratInputValidation:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    def strat_rast(self):
        return sgs.quantiles(self.rast, quantiles={'zq90': 3})

    def test_manual_weights_with_float_rounding(self):
        #np.sum([0.7, 0.2, 0.1]) == 0.9999999999999999
        samples = sgs.strat(self.strat_rast(), num_samples=30, allocation="manual", weights=[0.7, 0.2, 0.1], method="random")
        assert len(samples.samples_as_wkt()) > 0

    @pytest.mark.parametrize("method", ["random", "Queinnec"])
    def test_out_of_range_strata_value_raises(self, method):
        srast = self.strat_rast()
        srast.is_strat_rast = False #ignore metadata, so the num_strata below is used
        with pytest.raises(RuntimeError, match="num_strata"):
            sgs.strat(srast, num_samples=30, num_strata=2, method=method)

    def test_float_strat_raster_rejected(self):
        #sraster.tif stores whole-number strata as Float32, like R/terra outputs often do
        with pytest.raises(RuntimeError, match="not supported"):
            sgs.strat(sgs.SpatialRaster(sraster_geotiff_path), num_samples=30, num_strata=5)
