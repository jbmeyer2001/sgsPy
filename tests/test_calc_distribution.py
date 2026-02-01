import pytest
import numpy as np
import geopandas as gpd

import sgspy as sgs

from files import (
    mraster_geotiff_path,
    existing_shapefile_path 
)

class TestDistribution:
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    samples = sgs.SpatialVector(existing_shapefile_path)

    #we're testing with numpy (but not using numpy in the sgs package) because the sgspy distribution function
    #should be able to calculate the distribution on very large raster images which wouldn't
    #be able to effectively fit within memory in a numpy array
    def check(self, result, arr, bins, samples=False):
        [bin_vals, counts] = result["population"]

        farr = np.ndarray.flatten(arr)
        farr = farr[~np.isnan(farr)]
        [check_counts, check_bin_vals] = np.histogram(farr, bins=bins)

        assert np.array_equal(counts, check_counts)
        np.testing.assert_almost_equal(bin_vals, check_bin_vals, decimal=5)

        if samples:
            gdf = self.samples.to_geopandas()
            sample_vals = []

            for point in gdf['geometry']:
                x = int((point.x - self.rast.xmin) / self.rast.pixel_width)
                y = int((self.rast.ymax - point.y) / self.rast.pixel_height)
                sample_vals.append(arr[y, x])

            [bin_vals, counts] = result["sample"]
            [check_counts, check_bin_vals] = np.histogram(sample_vals, bins=bins, range=(check_bin_vals[0], check_bin_vals[bins]))

            assert np.array_equal(counts, check_counts)

    def test_bin_number(self):
        bands = {
            'zq90': self.rast.band('zq90'),
            'pzabove2': self.rast.band('pzabove2'),
            'zsd': self.rast.band('zsd')
        }
        
        for bins in [1, 10, 50, 100]:
            for band in ['zq90', 'pzabove2', 'zq90']:
                result = sgs.calculate.distribution(self.rast, band=band, bins=bins, plot=False)
                self.check(result, bands[band], bins)

    def test_sample_dist(self):
        bands = {
            'zq90': self.rast.band('zq90'),
            'pzabove2': self.rast.band('pzabove2'),
            'zsd': self.rast.band('zsd')
        }

        bins = 50

        for band in ['zq90', 'pzabove2', 'zsd']:
            result = sgs.calculate.distribution(self.rast, band=band, bins=bins, samples=self.samples, plot=False)
            self.check(result, bands[band], bins, True)
