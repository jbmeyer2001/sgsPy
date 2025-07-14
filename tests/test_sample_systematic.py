import geopandas as gpd
import pytest

import sgs

from files import mraster_geotiff_path

class TestSystematic:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    def check_samples(self, samples):
        for sample in samples:
            assert sample.x >= self.rast.xmin
            assert sample.x <= self.rast.xmax
            assert sample.y >= self.rast.ymin
            assert sample.y <= self.rast.ymax

    def test_points_in_correct_area(self):
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "square", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "triangle", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "hexagon", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "square", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "triangle", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "hexagon", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "square", "random"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "triangle", "random"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "hexagon", "random"))
        self.check_samples(samples)

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "square", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "triangle", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "square", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "triangle", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "square", "random"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "triangle", "random"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "random"))
        self.check_samples(samples)

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "square", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "triangle", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "hexagon", "centers"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "square", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "triangle", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "hexagon", "corners"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "square", "random"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "triangle", "random"))
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "hexagon", "random"))
        self.check_samples(samples)

    def test_inputs(self):
        sgs.systematic(self.rast, 100)
        sgs.systematic(self.rast, 100, "square", "centers")
        sgs.systematic(self.rast, 100, "triangle", "corners")
        sgs.systematic(self.rast, 100, "hexagon", "random")

        with pytest.raises(ValueError):
            sgs.systematic(self.rast, -1, "square", "centers")

        with pytest.raises(ValueError):
            sgs.systematic(self.rast, 0, "square", "centers")

        with pytest.raises(ValueError):
            sgs.systematic(self.rast, 100, "", "centers")

        with pytest.raises(ValueError):
            sgs.systematic(self.rast, 100, "square", "")

        with pytest.raises(ValueError):
            sgs.systematic(self.rast, 100, "square", "center")

        with pytest.raises(ValueError):
            sgs.systematic(self.rast, 100, "squares", "centers")

    def test_write(self, tmp_path):
        temp_dir = tmp_path / "test_output"
        temp_dir.mkdir()
        temp_file = temp_dir / "vect.shp"

        gs_samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "centers", filename=str(temp_file)))
        gs_file = gpd.read_file(temp_file)

        assert len(gs_samples.intersection(gs_file)) == len(gs_samples)
