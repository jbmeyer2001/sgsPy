import geopandas as gpd
import numpy as np
import pytest

import sgs

from files import (
        mraster_geotiff_path,
        access_shapefile_path
)

class TestBalanced:
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    srast = sgs.stratify.quantiles(rast, num_strata={'zq90': 5})
    access = sgs.SpatialVector(access_shapefile_path)
    num_data_pixels = 91195 #constant which is true for mraster.tif

    def check_samples(self, samples):
        for point in samples:
            assert point.x <= self.rast.xmax
            assert point.x >= self.rast.xmin
            assert point.y <= self.rast.ymax
            assert point.y >= self.rast.ymin

            x = (point.x - self.rast.xmin) / self.rast.pixel_width
            y = self.rast.height - ((point.y - self.rast.ymin) / self.rast.pixel_height)

            for band in self.rast.bands:
                assert not np.isnan(self.rast[band, int(y), int(x)])

    def check_samples_access(self, samples, buff_inner, buff_outer):
        gs_access = gpd.read_file(access_shapefile_path)
        if (buff_inner == 0):
            accessable = gs_access.buffer(buff_outer).union_all()
        else:
            accessable = gs_access.buffer(buff_outer).union_all().difference(gs_access.buffer(buff_inner).union_all())

        for sample in samples:
            assert accessable.contains(sample)

    def test_lcube(self):
        #just lcube
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lcube"
        ).samples_as_wkt())
        self.check_samples(samples)

        #lcube with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=500,
            algorithm="lcube",
            access=self.access,
            layer_name='access',
            buff_inner=30,
            buff_outer=300
        ).samples_as_wkt())
        self.check_samples(samples)
        self.check_samples_access(samples, 30, 300)

        #lcube with prob
        prob = np.full(
            self.num_data_pixels,
            50/self.num_data_pixels,
            np.float64
        )
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lcube",
            prob=prob,
        ).samples_as_wkt())
        self.check_samples(samples)

        #lcube with prob not all equal
        prob = np.concatenate(
            (np.full(
                91195 - 2000,
                0,
                np.float64
            ),
            np.full(
                2000,
                50 / 2000,
                np.float64
            ))
        )
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lcube",
            prob=prob,
        ).samples_as_wkt())
        self.check_samples(samples)

    def test_lcubestratified(self):
        #just lcubestratified
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lcubestratified",
            srast=self.srast,
            srast_band='strat_zq90'
        ).samples_as_wkt())
        self.check_samples(samples)

        #lcubestratified with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=500,
            algorithm="lcubestratified",
            srast=self.srast,
            srast_band='strat_zq90',
            access=self.access,
            layer_name='access',
            buff_inner=30,
            buff_outer=300
        ).samples_as_wkt())
        self.check_samples(samples)
        self.check_samples_access(samples, 30, 300)

        #lcubestratified with prob
        prob = np.full(
            self.num_data_pixels,
            50/self.num_data_pixels,
            np.float64
        )
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lcubestratified",
            srast=self.srast,
            srast_band='strat_zq90',
            prob=prob,
        ).samples_as_wkt())
        self.check_samples(samples)

        #lcubestratified with prob not all equal
        prob = np.concatenate(
            (np.full(
                91195 - 2000,
                0,
                np.float64
            ),
            np.full(
                2000,
                50 / 2000,
                np.float64
            ))
        )
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lcubestratified",
            srast=self.srast,
            srast_band='strat_zq90',
            prob=prob,
        ).samples_as_wkt())
        self.check_samples(samples)

    def test_lpm2_kdtree(self):
        #just lpm2_kdtree
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lpm2_kdtree"
        ).samples_as_wkt())
        self.check_samples(samples)

        #lpm2_kdtree with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=500,
            algorithm="lpm2_kdtree",
            access=self.access,
            layer_name='access',
            buff_inner=30,
            buff_outer=300
        ).samples_as_wkt())
        self.check_samples(samples)
        self.check_samples_access(samples, 30, 300)

        #lpm2_kdtree with prob
        prob = np.full(
            self.num_data_pixels,
            50/self.num_data_pixels,
            np.float64
        )
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lpm2_kdtree",
            prob=prob,
        ).samples_as_wkt())
        self.check_samples(samples)

        #lpm2_kdtree with prob not all equal
        prob = np.concatenate(
            (np.full(
                91195 - 2000,
                0,
                np.float64
            ),
            np.full(
                2000,
                50 / 2000,
                np.float64
            ))
        )
        samples = gpd.GeoSeries.from_wkt(sgs.sample.balanced(
            self.rast,
            num_samples=1000,
            algorithm="lpm2_kdtree",
            prob=prob,
        ).samples_as_wkt())
        self.check_samples(samples)

    def test_inputs(self):
        for bad_algorithm in ["lpm_kdtree", "lcube_stratified", "", "cube"]:
            with pytest.raises(ValueError):
                samples = sgs.sample.balanced(self.rast, num_samples=500, algorithm=bad_algorithm)

        for bad_samples in [0, -1000]:
            with pytest.raises(ValueError):
                samples = sgs.sample.balanced(self.rast, num_samples=bad_samples)

        with pytest.raises(ValueError):
            samples = sgs.sample.balanced(self.rast, num_samples=500, algorithm="lcubestratified")

        with pytest.raises(ValueError):
            samples = sgs.sample.balanced(self.rast, num_samples=500, algorithm="lcubestratified", srast=self.srast)

        #TODO access checks??
