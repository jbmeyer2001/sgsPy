import geopandas as gpd
import numpy as np
import pytest

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    mraster_small_geotiff_path,
    access_shapefile_path,
    existing_shapefile_path,
)

class TestSrs:
    rast = sgs.SpatialRaster(mraster_small_geotiff_path)
    
    mrast_full = sgs.SpatialRaster(mraster_geotiff_path)
    access = sgs.SpatialVector(access_shapefile_path)
    existing = sgs.SpatialVector(existing_shapefile_path)

    def test_num_points(self):
        with pytest.raises(ValueError):
            samples = sgs.srs(self.rast, num_samples=0)

        sample = sgs.srs(self.rast, num_samples=1).samples_as_wkt()
        assert len(sample) == 1

        samples = sgs.srs(self.rast, num_samples=50).samples_as_wkt()
        assert len(samples) == 50

        samples = sgs.srs(self.rast, num_samples=2000).samples_as_wkt()
        assert len(samples) == 2000

    def test_mindist(self):
        def check_samples(mindist, samples):
            gs = gpd.GeoSeries.from_wkt(samples)
            for i in range(len(samples) - 1):
                distances = gs[i].distance(gs[i+1:])
                assert 0 == int(np.sum(np.array([distances < mindist]).astype(int)))

        mindist = 1
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000).samples_as_wkt() 
        check_samples(mindist, samples)

        mindist = 50
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000).samples_as_wkt()
        check_samples(mindist, samples)

        mindist = 200.5
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000).samples_as_wkt()
        check_samples(mindist, samples)

        mindist = 6000
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000).samples_as_wkt()
        check_samples(mindist, samples)

    def test_points_in_bounds(self):
        samples = sgs.srs(self.rast, num_samples=1000).samples_as_wkt()
        gs = gpd.GeoSeries.from_wkt(samples)
        
        #ensure all points are within raster bounds
        for point in gs:
           assert point.x <= self.rast.xmax
           assert point.x >= self.rast.xmin
           assert point.y <= self.rast.ymax
           assert point.y >= self.rast.ymin

    def test_points_not_nan(self):
        samples = sgs.srs(self.rast, num_samples=1000).samples_as_wkt()
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
 
    def test_write_output(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "vect.shp"
        gs_samples = sgs.srs(self.rast, num_samples = 1000, filename=str(temp_file)).to_geopandas()['geometry']
        gs_file = gpd.read_file(temp_file)
        
        assert len(gs_samples.intersection(gs_file)) == 1000

    def test_access(self):
        gs_access = gpd.read_file(access_shapefile_path)

        #test just buff_outer working
        samples = gpd.GeoSeries.from_wkt(sgs.srs(self.mrast_full, 50000, access=self.access, buff_outer=100).samples_as_wkt())
        accessable = gs_access.buffer(100).union_all()
        for sample in samples:
            assert accessable.contains(sample)

        #test buff_outer works with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.srs(self.mrast_full, 50000, 200, access=self.access, buff_outer=100).samples_as_wkt())
        #accessable stays the same because buff_outer is the same
        for sample in samples:
            assert accessable.contains(sample)

        #test buff_outer and buff_inner works
        bad_points=[]
        samples = gpd.GeoSeries.from_wkt(sgs.srs(self.mrast_full, 50000, access=self.access, buff_outer=200, buff_inner=100).samples_as_wkt())
        accessable = gs_access.buffer(200).union_all().difference(gs_access.buffer(100).union_all())
        for sample in samples:
            assert accessable.contains(sample)

        samples = gpd.GeoSeries.from_wkt(sgs.srs(self.mrast_full, 50000, 200, access=self.access, buff_outer=200, buff_inner=100).samples_as_wkt())
        #accessable stays the same because buff_outer and buff_inner are the same
        for sample in samples:
            assert accessable.contains(sample)

    def test_existing(self):
        existing = gpd.read_file(existing_shapefile_path)['geometry']
        samples = gpd.GeoSeries.from_wkt(sgs.srs(self.mrast_full, 1000, existing=self.existing).samples_as_wkt())

        for sample in existing:
            assert samples.contains(sample).any()
    
    #TODO test input values
