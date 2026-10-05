import time

import geopandas as gpd
import shapely

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    existing_shapefile_path,
    access_shapefile_path,
)


class TestSystematicEdgeCases:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    def test_corners_with_existing_has_no_duplicates(self):
        """
        systematic.h: with existing, the second corner used to be written as a copy of the first one.
        Hexagon corners never coincide, so a duplicate within one run means a corner was written twice.
        """
        existing = sgs.SpatialVector(existing_shapefile_path)
        samples = gpd.GeoSeries.from_wkt(sgs.systematic(self.rast, 300, "hexagon", "corners", existing=existing).samples_as_wkt())
        assert not samples.duplicated().any()
