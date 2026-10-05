import sys

import numpy as np
import pytest

import sgspy as sgs

from files import mraster_small_geotiff_path

#mraster_small.tif is 141 x 110 in 4-row strips: 28 block rows, the last one holding only 2 rows
MANY_THREADS = 32


def force_large_raster(monkeypatch, module):
    monkeypatch.setattr(sys.modules["sgspy." + module], "GIGABYTE", 1)


def small_rast():
    return sgs.SpatialRaster(mraster_small_geotiff_path)


class TestLargeRasterMatchesInMemory:
    def test_pca(self, monkeypatch):
        expected = sgs.pca(small_rast(), num_comp=2)
        force_large_raster(monkeypatch, "calculate.pca.pca")
        result = sgs.pca(small_rast(), num_comp=2)

        for i in range(2):
            a, b = expected.band(i), result.band(i)
            assert np.array_equal(np.isnan(a), np.isnan(b))
            #eigenvectors are only defined up to sign
            np.testing.assert_allclose(np.abs(a), np.abs(b), atol=1e-4, equal_nan=True)

    @pytest.mark.parametrize("map", [False, True])
    def test_breaks(self, monkeypatch, map):
        kw = dict(breaks={'zq90': [3, 5, 11, 18], 'zsd': [2, 5]}, map=map)
        expected = sgs.breaks(small_rast(), **kw)
        force_large_raster(monkeypatch, "stratify.breaks.breaks")
        result = sgs.breaks(small_rast(), thread_count=MANY_THREADS, **kw)

        for band in expected.bands:
            assert np.array_equal(expected.band(band), result.band(band)), band

    @pytest.mark.parametrize("map", [False, True])
    def test_quantiles(self, monkeypatch, map):
        kw = dict(quantiles={'zq90': 4, 'zsd': 4}, map=map)
        expected = sgs.quantiles(small_rast(), **kw)
        force_large_raster(monkeypatch, "stratify.quantiles.quantiles")
        result = sgs.quantiles(small_rast(), thread_count=MANY_THREADS, **kw)

        #the large path uses MKL streaming quantiles, which are approximate
        #so we give ourselves a bit of slack here
        for band in expected.bands:
            assert np.mean(expected.band(band) != result.band(band)) < 0.005, band

    def test_map(self, monkeypatch):
        strat = sgs.breaks(small_rast(), breaks={'zq90': [3, 5, 11, 18], 'zsd': [2, 5]})
        args = (strat, ['strat_zq90', 'strat_zsd'], [5, 3])
        expected = sgs.map(args)
        force_large_raster(monkeypatch, "stratify.map.map")
        result = sgs.map(args, thread_count=MANY_THREADS)

        assert np.array_equal(expected.band(0), result.band(0))
