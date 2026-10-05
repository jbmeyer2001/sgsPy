import geopandas as gpd
import numpy as np

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    existing_shapefile_path,
    access_shapefile_path,
)


def points(xs, ys):
    crs = gpd.read_file(access_shapefile_path).crs
    return sgs.SpatialVector.from_geopandas(gpd.GeoDataFrame(geometry=gpd.points_from_xy(xs, ys), crs=crs))


class TestClhsExisting:
    rast = sgs.SpatialRaster(mraster_geotiff_path)

    def test_new_samples_avoid_bins_covered_by_existing(self):
        num_samples, num_existing = 40, 20
        zq90 = self.rast.band('zq90')
        valid = np.all([~np.isnan(self.rast.band(b)) for b in self.rast.bands], axis=0)
        edges = np.nanquantile(zq90[valid], np.linspace(0, 1, num_samples + 1))

        #one existing plot in each of the top 20 (of 40) zq90 quantile bins
        rows, cols = [], []
        for b in range(num_samples - num_existing, num_samples):
            r, c = np.nonzero(valid & (zq90 > edges[b]) & (zq90 < edges[b + 1]))
            rows.append(r[len(r) // 2])
            cols.append(c[len(c) // 2])
        rows, cols = np.array(rows), np.array(cols)
        xs = self.rast.xmin + (cols + 0.5) * self.rast.pixel_width
        ys = self.rast.ymax - (rows + 0.5) * self.rast.pixel_height
        crs = gpd.read_file(access_shapefile_path).crs
        existing = sgs.SpatialVector.from_geopandas(gpd.GeoDataFrame(geometry=gpd.points_from_xy(xs, ys), crs=crs))

        fractions = []
        for _ in range(3):
            pts = gpd.GeoSeries.from_wkt(sgs.clhs(self.rast, num_samples=num_samples, existing=existing).samples_as_wkt())
            pr = ((self.rast.ymax - pts.y) / self.rast.pixel_height).astype(int)
            pc = ((pts.x - self.rast.xmin) / self.rast.pixel_width).astype(int)
            is_existing = np.isin(pr * zq90.shape[1] + pc, rows * zq90.shape[1] + cols)
            assert is_existing.sum() == num_existing #replace=0 keeps every existing plot
            new_values = zq90[pr[~is_existing], pc[~is_existing]]
            fractions.append(np.mean(new_values > edges[num_samples - num_existing]))

        #ignoring the existing plots puts ~half of the new samples into the covered bins
        assert np.mean(fractions) < 0.25

    def test_more_existing_than_num_samples(self):
        existing = sgs.SpatialVector(existing_shapefile_path) #200 points
        samples = sgs.clhs(self.rast, num_samples=50, existing=existing)
        assert len(samples.samples_as_wkt()) == 200

    def test_existing_on_nan_or_shared_pixel_is_dropped(self):
        #existing.h keeps one point per pixel; readRaster only adds it if the pixel is valid
        valid = np.all([~np.isnan(self.rast.band(b)) for b in self.rast.bands], axis=0)
        (nr, nc), (vr, vc) = np.argwhere(~valid)[0], np.argwhere(valid)[0]
        off = np.array([.5, .5, .25])
        xs = self.rast.xmin + (np.array([nc, vc, vc]) + off) * self.rast.pixel_width
        ys = self.rast.ymax - (np.array([nr, vr, vr]) + off) * self.rast.pixel_height
        samples = sgs.clhs(self.rast, 10, existing=points(xs, ys)).to_geopandas()
        assert (samples['existing'] == 1).sum() == 1
        assert len(samples) == 10
