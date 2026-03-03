import warnings

import geopandas as gpd
import numpy as np
import pytest

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    access_shapefile_path,
    existing_shapefile_path,
)

class TestClhs:
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    access = sgs.SpatialVector(access_shapefile_path)
    existing = sgs.SpatialVector(existing_shapefile_path)

    def test_num_points(self):
        with pytest.raises(ValueError):
            samples = sgs.sample.clhs(self.rast, num_samples=0)

        sample = sgs.sample.clhs(self.rast, num_samples=10).samples_as_wkt()
        assert len(sample) == 10

        sample = sgs.sample.clhs(self.rast, num_samples=100).samples_as_wkt()
        assert len(sample) == 100

        sample = sgs.sample.clhs(self.rast, num_samples=250).samples_as_wkt()
        assert len(sample) == 250

    def test_points_in_bounds(self):
        for _ in range(10):
            samples = sgs.sample.clhs(self.rast, num_samples=200).samples_as_wkt()
            gs = gpd.GeoSeries.from_wkt(samples)

            for point in gs:
                assert point.x <= self.rast.xmax
                assert point.x >= self.rast.xmin
                assert point.y <= self.rast.ymax
                assert point.y >= self.rast.ymin

    def test_points_not_nan(self):
        for _ in range(10):
            samples = sgs.sample.clhs(self.rast, num_samples=200).samples_as_wkt()
            gs = gpd.GeoSeries.from_wkt(samples)

            #find indexes for all points and ensure they're not nan pixels.
            for point in gs:
                x_index = (point.x - self.rast.xmin) / self.rast.pixel_width
                y_index = self.rast.height - ((point.y - self.rast.ymin) / self.rast.pixel_height) #origin at top left instead of bottom left
                pixel_value = self.rast.band(0)[int(y_index), int(x_index)]
                if np.isnan(pixel_value):
                    print(point.x)
                    print(point.y)
                    print(x_index)
                    print(y_index)
                    assert False

    def test_access(self):
        gs_access = gpd.read_file(access_shapefile_path)
       
        #both access tests, one with just buff_outer and one with both buff_outer and buff_inner
        accessible1 = gs_access.buffer(100).union_all()
        accessible2 = gs_access.buffer(200).union_all().difference(gs_access.buffer(100).union_all())

        for _ in range(5):
            #just buff_outer
            samples = gpd.GeoSeries.from_wkt(
                sgs.sample.clhs(self.rast, 200, access=self.access, buff_outer=100).samples_as_wkt()
            )
            for sample in samples:
                assert accessible1.contains(sample)

            #both buff_inner and buff_outer
            samples = gpd.GeoSeries.from_wkt(
                sgs.sample.clhs(self.rast, 200, access=self.access, buff_outer=200, buff_inner=100).samples_as_wkt()
            )
            for sample in samples:
                assert accessible2.contains(sample)

    def test_existing(self):
        existing_set = set(gpd.read_file(existing_shapefile_path)['geometry'])

        #test replace negative
        with pytest.raises(ValueError):
            samples = sgs.sample.clhs(self.rast, 200, existing=self.existing, replace = -1)

        #test replace 0
        for replace in [0, 10, 50, 100]:
            samples = sgs.sample.clhs(self.rast, 200, existing=self.existing, replace=replace).to_geopandas()

            replaced = samples[samples["existing"] == 0]
            assert len(replaced['geometry']) <= replace

            kept = set(samples[samples["existing"] == 1]['geometry'])
            assert len(kept.difference(existing_set)) == 0

        #test 50 required samples after existing
        samples = sgs.sample.clhs(self.rast, 250, existing=self.existing).to_geopandas()
        samples_set = set(samples['geometry'])

        assert len(existing_set.difference(samples_set)) == 0
        assert len(samples_set.difference(existing_set)) == 50

    def test_write(self, tmp_path):
        temp_dir = tmp_path / "test_output"
        temp_dir.mkdir()
        temp_file = temp_dir / "vect.shp"

        gs_samples = sgs.sample.clhs(self.rast, 200, filename=str(temp_file)).to_geopandas()['geometry']
        gs_file = gpd.read_file(temp_file)

        # on linux, the following throws a warning about the spatial references differing by:
        # UTM Zone 17, Northern Hemisphere vs UTM_Zone_17_Northern_Hemisphere
        # which shouldn't cause any kind of problem
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert len(gs_samples.intersection(gs_file)) == len(gs_samples)

