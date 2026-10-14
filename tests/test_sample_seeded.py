import sgspy as sgs
import pytest

from files import (
    mraster_geotiff_path,
)


class TestSeeded:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    def test_negative_seed(self):
        with pytest.raises(ValueError, match="Expected non-negative integer"):
            sgs.srs(self.rast, num_samples=50, random_state=-1)

    def test_random_state_is_reproducible_srs(self):
        run = lambda seed: sorted(sgs.srs(self.rast, num_samples=50, random_state=seed).samples_as_wkt())
        assert run(7) == run(7)
        assert run(7) != run(8)
        assert run(None) != run(None)

    def test_random_state_is_reproducible_strat(self):
        """Stratify otherwise Float32 is SEGFAULT"""
        srast = sgs.stratify.quantiles(self.rast, quantiles={"zq90": 6})
        run = lambda seed: sorted(sgs.strat(srast, num_strata=6, num_samples=50, random_state=seed).samples_as_wkt())
        assert run(9) == run(9)
        assert run(9) != run(10)
        assert run(None) != run(None)

    def test_random_state_is_reproducible_sys(self):
        run = lambda seed: sorted(sgs.systematic(self.rast, cellsize=100, random_state=seed).samples_as_wkt())
        assert run(11) == run(11)
        assert run(11) != run(12)
        assert run(None) != run(None)

    def test_random_state_is_reproducible_clhs(self):
        run = lambda seed: sorted(sgs.clhs(self.rast, num_samples=50, random_state=seed).samples_as_wkt())
        assert run(13) == run(13)
        assert run(13) != run(14)
        assert run(None) != run(None)
