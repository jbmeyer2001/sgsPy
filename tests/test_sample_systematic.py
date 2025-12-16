import numpy as np
import geopandas as gpd
import pytest

import sgs

from files import (
    mraster_geotiff_path,
    access_shapefile_path,
    existing_shapefile_path
)

class TestSystematic:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    access = sgs.SpatialVector(access_shapefile_path)
    existing = sgs.SpatialVector(existing_shapefile_path)

    def check_samples(self, samples):
        for sample in samples:
            assert sample.x >= self.rast.xmin
            assert sample.x <= self.rast.xmax
            assert sample.y >= self.rast.ymin
            assert sample.y <= self.rast.ymax

    def test_points_in_correct_area(self):
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "square", "centers").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "hexagon", "centers").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "square", "corners").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "hexagon", "corners").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "square", "random").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 30, "hexagon", "random").samples_as_wkt())
        self.check_samples(samples)

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "square", "centers").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "centers").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "square", "corners").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "corners").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "square", "random").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 500, "hexagon", "random").samples_as_wkt())
        self.check_samples(samples)

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "square", "centers").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "hexagon", "centers").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "square", "corners").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "hexagon", "corners").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "square", "random").samples_as_wkt())
        self.check_samples(samples)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 3000, "hexagon", "random").samples_as_wkt())
        self.check_samples(samples)

    def test_inputs(self):
        sgs.systematic(self.rast, 100)
        sgs.systematic(self.rast, 100, "square", "centers")
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

        gs_samples = sgs.systematic(self.rast, 500, "hexagon", "centers", filename=str(temp_file)).to_geopandas()['geometry']
        gs_file = gpd.read_file(temp_file)

        assert len(gs_samples.intersection(gs_file)) == len(gs_samples)

    def test_existing(self):
        existing = gpd.read_file(existing_shapefile_path)['geometry']

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "square", "centers", existing=self.existing).samples_as_wkt())
        for sample in existing: assert samples.contains(sample).any()

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "centers", existing=self.existing).samples_as_wkt())
        for sample in existing: assert samples.contains(sample).any()
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "square", "corners", existing=self.existing).samples_as_wkt())
        for sample in existing: assert samples.contains(sample).any()
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "corners", existing=self.existing).samples_as_wkt())
        for sample in existing: assert samples.contains(sample).any()
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "square", "random", existing=self.existing).samples_as_wkt())
        for sample in existing: assert samples.contains(sample).any()
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "random", existing=self.existing).samples_as_wkt())
        for sample in existing: assert samples.contains(sample).any()


    def test_access(self):
        gs_access = gpd.read_file(access_shapefile_path)
        accessible = gs_access.buffer(120).union_all()

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "square", "centers", access=self.access, buff_outer=120).samples_as_wkt())
        for sample in samples: assert accessible.contains(sample)

        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "centers", access=self.access, buff_outer=120).samples_as_wkt())
        for sample in samples: assert accessible.contains(sample)
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "square", "corners", access=self.access, buff_outer=120).samples_as_wkt())
        for sample in samples: assert accessible.contains(sample)
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "corners", access=self.access, buff_outer=120).samples_as_wkt())
        for sample in samples: assert accessible.contains(sample)
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "square", "random", access=self.access, buff_outer=120).samples_as_wkt())
        for sample in samples: assert accessible.contains(sample)
        
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "random", access=self.access, buff_outer=120).samples_as_wkt())
        for sample in samples: assert accessible.contains(sample)

