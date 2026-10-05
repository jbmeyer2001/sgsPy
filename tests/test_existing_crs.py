import geopandas as gpd
import pytest

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    existing_shapefile_path,
)


class TestExistingCrsMismatch:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    def wrong_crs_existing(self):
        gdf = gpd.read_file(existing_shapefile_path).set_crs("EPSG:4326", allow_override=True)
        return sgs.SpatialVector.from_geopandas(gdf)

    @pytest.mark.parametrize("sampler", [
        lambda rast, existing: sgs.srs(rast, num_samples=10, existing=existing),
        lambda rast, existing: sgs.clhs(rast, num_samples=10, existing=existing),
        lambda rast, existing: sgs.systematic(rast, 500, "square", "centers", existing=existing),
    ], ids=["srs", "clhs", "systematic"])
    def test_raises(self, sampler):
        with pytest.raises(RuntimeError, match="spatial reference"):
            sampler(self.rast, self.wrong_crs_existing())
