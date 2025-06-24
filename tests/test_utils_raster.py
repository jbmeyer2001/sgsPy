import pytest
import sgs
import numpy as np
import matplotlib.pyplot as plt

from files import (
    mraster_geotiff_path,
    mraster_small_geotiff_path,
    sraster_geotiff_path,
    sraster2_geotiff_path,
    mraster_zq90_path,
    mraster_pzabove2_path,
    mraster_zsd_path,
    mraster_small_zq90_path,
    mraster_small_pzabove2_path,
    mraster_small_zsd_path,
    sraster_strata_path,
    sraster2_band_path,
)

class TestSpatialRaster:
    def mraster_check(self, rast):
        assert rast.width == 373
        assert rast.height == 277
        assert rast.band_count == 3
        assert rast.pixel_width == 20
        assert rast.pixel_height == 20
        assert rast.xmin == 431100
        assert rast.xmax == 438560
        assert rast.ymin == 5337700
        assert rast.ymax == 5343240
        assert 'zq90' in rast.bands
        assert 'pzabove2' in rast.bands
        assert 'zsd' in rast.bands
        assert len(rast.bands) == 3
        assert np.array_equal(np.load(mraster_zq90_path), rast[0], equal_nan=True)
        assert np.array_equal(np.load(mraster_pzabove2_path), rast[1], equal_nan=True)
        assert np.array_equal(np.load(mraster_zsd_path), rast[2], equal_nan=True)

    def mraster_small_check(self, rast):
        assert rast.width == 141
        assert rast.height == 110
        assert rast.band_count == 3
        assert rast.pixel_width == 20
        assert rast.pixel_height == 20
        assert rast.xmin == 432440
        assert rast.xmax == 435260
        assert rast.ymin == 5340000
        assert rast.ymax == 5342200
        assert 'zq90' in rast.bands
        assert 'pzabove2' in rast.bands
        assert 'zsd' in rast.bands
        assert len(rast.bands) == 3
        assert np.array_equal(np.load(mraster_small_zq90_path), rast[0], equal_nan=True)
        assert np.array_equal(np.load(mraster_small_pzabove2_path), rast[1], equal_nan=True)
        assert np.array_equal(np.load(mraster_small_zsd_path), rast[2], equal_nan=True)
    
    def sraster_check(self, rast):
        assert rast.width == 373
        assert rast.height == 277
        assert rast.band_count == 1
        assert rast.pixel_width == 20
        assert rast.pixel_height == 20
        assert rast.xmin == 431100
        assert rast.xmax == 438560
        assert rast.ymin == 5337700
        assert rast.ymax == 5343240
        assert 'strata' in rast.bands
        assert len(rast.bands) == 1
        assert np.array_equal(np.load(sraster_strata_path), rast[0], equal_nan=True)

    def sraster2_check(self, rast):
        assert rast.width == 861
        assert rast.height == 611
        assert rast.band_count == 1
        assert rast.pixel_width - 0.0003353943 == pytest.approx(0, abs=1e-10)
        assert rast.pixel_height - 0.0003353943 == pytest.approx(0, abs=1e-10)
        assert rast.xmin - -71.9365 == pytest.approx(0, abs=1e-4)
        assert rast.xmax - -71.64773 == pytest.approx(0, abs=1e-4)
        assert rast.ymin - 45.4883 == pytest.approx(0, abs=1e-4)
        assert rast.ymax - 45.69332 == pytest.approx(0, abs=1e-4)
        assert '' in rast.bands
        assert len(rast.bands) == 1
        assert np.array_equal(np.load(sraster2_band_path), rast[0], equal_nan=True)

    def test_construct_from_path(self):
        rast = sgs.utils.raster.SpatialRaster(mraster_geotiff_path)
        self.mraster_check(rast)

        rast = sgs.utils.raster.SpatialRaster(mraster_small_geotiff_path)
        self.mraster_small_check(rast)

        rast = sgs.utils.raster.SpatialRaster(sraster_geotiff_path)
        self.sraster_check(rast)

        rast = sgs.utils.raster.SpatialRaster(sraster2_geotiff_path)
        self.sraster2_check(rast)

    def test_construct_from_existing(self):
        rast = sgs.utils.raster.SpatialRaster(mraster_geotiff_path)
        new_rast = sgs.utils.raster.SpatialRaster(rast.cpp_raster)
        self.mraster_check(new_rast)

        rast = sgs.utils.raster.SpatialRaster(mraster_small_geotiff_path)
        new_rast = sgs.utils.raster.SpatialRaster(rast.cpp_raster)
        self.mraster_small_check(new_rast)

        rast = sgs.utils.raster.SpatialRaster(sraster_geotiff_path)
        new_rast = sgs.utils.raster.SpatialRaster(rast.cpp_raster)
        self.sraster_check(new_rast)

        rast = sgs.utils.raster.SpatialRaster(sraster2_geotiff_path)
        new_rast = sgs.utils.raster.SpatialRaster(rast.cpp_raster)
        self.sraster2_check(new_rast)

    def test_raster_slicing(self):
        rast = sgs.utils.raster.SpatialRaster(mraster_small_geotiff_path)
        zq90 = np.load(mraster_small_zq90_path)
        pzabove2 = np.load(mraster_small_pzabove2_path)
        zsd = np.load(mraster_small_zsd_path)

        assert np.array_equal(zq90, rast[0], equal_nan=True)
        assert np.array_equal(pzabove2, rast[1], equal_nan=True)
        assert np.array_equal(zsd, rast[2], equal_nan=True)

        assert np.array_equal(zq90, rast['zq90'], equal_nan=True)
        assert np.array_equal(pzabove2, rast['pzabove2'], equal_nan=True)
        assert np.array_equal(zsd, rast['zsd'], equal_nan=True)

        assert np.array_equal(zq90[0:10, 0:10], rast['zq90', 0:10, 0:10], equal_nan=True)
        assert np.array_equal(zq90[131:141, 100:110], rast['zq90', 131:141, 100:110], equal_nan=True)
        assert np.array_equal([zq90, pzabove2], rast[0:2], equal_nan=True)
        assert np.array_equal(zq90[73, 46], rast['zq90', 73, 46], equal_nan=True)
