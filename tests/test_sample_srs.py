import geopandas as gpd
import numpy as np
import pytest

import sgs

from files import mraster_small_geotiff_path

class TestSrs:
    rast = sgs.SpatialRaster(mraster_small_geotiff_path)

    def test_num_points(self):
        with pytest.raises(ValueError):
            samples = sgs.srs(self.rast, num_samples=0)

        samples = sgs.srs(self.rast, num_samples=1)
        assert len(samples) == 1

        samples = sgs.srs(self.rast, num_samples=50)
        assert len(samples) == 50

        samples = sgs.srs(self.rast, num_samples=2000)
        assert len(samples) == 2000

        with pytest.raises(RuntimeError):
            samples = sgs.srs(self.rast, num_samples=200000)

    def test_mindist(self):
        def check_samples(mindist, samples):
            gs = gpd.GeoSeries.from_wkt(samples)
            for i in range(len(samples) - 1):
                distances = gs[i].distance(gs[i+1:])
                assert 0 == int(np.sum(np.array([distances < mindist]).astype(int)))

        mindist = 1
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000) 
        check_samples(mindist, samples)

        mindist = 50
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000) 
        check_samples(mindist, samples)

        mindist = 200.5
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000) 
        check_samples(mindist, samples)

        mindist = 6000
        samples = sgs.srs(self.rast, mindist=mindist, num_samples=1000) 
        check_samples(mindist, samples)

    def test_points_in_bounds(self):
        samples = sgs.srs(self.rast, num_samples=1000)
        gs = gpd.GeoSeries.from_wkt(samples)
        
        #ensure no points selected cover nan pixels
        for point in gs:
           assert point.x <= self.rast.xmax
           assert point.x >= self.rast.xmin
           assert point.y <= self.rast.ymax
           assert point.y >= self.rast.ymin

    def test_points_not_nan(self):
        samples = sgs.srs(self.rast, num_samples=1000)
        gs = gpd.GeoSeries.from_wkt(samples)

        #find indexes for all points and ensure they're not nan pixels.
        for point in gs:
            x_index = (point.x - self.rast.xmin) / self.rast.pixel_width
            y_index = self.rast.height - ((point.y - self.rast.ymin) / self.rast.pixel_height) #origin at top left instead of bottom left
            pixel_value = self.rast[0, int(y_index), int(x_index)]
            if np.isnan(pixel_value):
                print(point.x)
                print(point.y)
                print(x_index)
                print(y_index)
                assert False
 
    def test_output(self, tmp_path):
        temp_dir = tmp_path / "test_out"
        temp_dir.mkdir()

        temp_file = temp_dir / "vect.shp"
        samples = sgs.srs(self.rast, num_samples = 1000, filename=str(temp_file))
        gs_samples = gpd.GeoSeries.from_wkt(samples)
        gs_file = gpd.read_file(temp_file)
        
        assert len(gs_samples.intersection(gs_file)) == 1000

