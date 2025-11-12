import platform
import geopandas as gpd
import numpy as np
import pytest

import sgs

from files import (
    mraster_geotiff_path,
    access_shapefile_path,
    existing_shapefile_path,
)

class TestStrat:
    rast = sgs.SpatialRaster(mraster_geotiff_path)
    access = sgs.SpatialVector(access_shapefile_path)
    existing = sgs.SpatialVector(existing_shapefile_path)

    def check_points_in_bounds(self, srast, samples):
        for point in samples:
            assert point.x <= self.rast.xmax
            assert point.x >= self.rast.xmin
            assert point.y <= self.rast.ymax
            assert point.y >= self.rast.ymin

            #also check to ensure pixel isn't nan
            x = (point.x - self.rast.xmin) / self.rast.pixel_width
            y = self.rast.height - ((point.y - self.rast.ymin) / self.rast.pixel_height)

            pixel_val = srast.band(0)[int(y), int(x)]

            if np.isnan(pixel_val):
                assert False

    def check_mindist(self, mindist, samples):
        for i in range(len(samples) - 1):
            distances = samples[i].distance(samples[i+1:])
            assert 0 == int(np.sum(np.array([distances < mindist]).astype(int)))

    def check_access(self, samples, buff_inner, buff_outer):
        gs_access = gpd.read_file(access_shapefile_path)
        if (buff_inner == 0):
            accessable = gs_access.buffer(buff_outer).union_all()
        else:
            accessable = gs_access.buffer(buff_outer).union_all().difference(gs_access.buffer(buff_inner).union_all())

        for sample in samples:
            assert accessable.contains(sample)

    def check_focal_window(self, srast, samples, wrow, wcol):
        for point in samples:
            x = (point.x - self.rast.xmin) / self.rast.pixel_width
            y = self.rast.height - ((point.y - self.rast.ymin) / self.rast.pixel_height) #origin at top left

            for row in range(wrow):
                for col in range(wcol):
                    y_check = y + row - wrow // 2
                    x_check = x + col - wcol // 2
                    assert(srast.band(0)[int(y), int(x)] == srast.band(0)[int(y_check), int(x_check)])

    def get_allocation_percentages(self, srast, samples):
        allocation = {}
        for point in samples:
            x = (point.x - self.rast.xmin) / self.rast.pixel_width
            y = self.rast.height - ((point.y - self.rast.ymin) / self.rast.pixel_height)
            strata = srast.band(0)[int(y), int(x)]
           
            if strata in allocation:
                allocation[strata] += 1
            else:
                allocation[strata] = 1

        for key in allocation.keys():
            allocation[key] = allocation[key] / len(samples)

        return allocation

    def test_random_allocation_equal(self):
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 5})

        #without mindist or access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast, 
            'strat_zq90',
            num_samples=500,
            num_strata=5,
            allocation="equal",
            method="random",
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0)

        #with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=5,
            allocation="equal",
            method="random",
            mindist=150,
        ).samples_as_wkt())
        
        assert len(samples) > 490 #mindist means we might not get the full 500
        self.check_points_in_bounds(srast, samples)
        self.check_mindist(150, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0, abs=0.03)

        #with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=5,
            allocation="equal",
            method = "random",
            access = self.access,
            layer_name = 'access',
            buff_inner = 90,
            buff_outer = 300,
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 90, 300)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0, abs=0.03)

        #with access and mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=5,
            allocation="equal",
            method="random",
            access = self.access,
            layer_name = 'access',
            buff_outer=600,
            mindist=90,
        ).samples_as_wkt())


        assert len(samples) > 90 #mindist means we may not get the full 100
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 0, 600)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0, abs=0.03)

    def test_random_allocation_proportional(self):
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 8})

        #without mindist or access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast, 
            'strat_zq90',
            num_samples=500,
            num_strata=8,
            allocation="prop",
            method="random"
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

        #with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=8,
            allocation="prop",
            method="random",
            mindist=150,
        ).samples_as_wkt())
        
        assert len(samples) > 490 #mindist means we may not get the full 500
        self.check_points_in_bounds(srast, samples)
        self.check_mindist(150, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

        #with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=8,
            allocation="prop",
            method = "random",
            access = self.access,
            layer_name = 'access',
            buff_inner = 90,
            buff_outer = 300,
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 90, 300)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

        #with access and mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=128,
            num_strata=8,
            allocation="prop",
            method="random",
            access = self.access,
            layer_name = 'access',
            buff_outer=600,
            mindist=90,
        ).samples_as_wkt())

        assert len(samples) > 115 #mindist means we may not get the full 128 
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 0, 600)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

    def test_random_allocation_manual(self):
        #without mindist or access
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 4})
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=4,
            allocation="manual",
            method="random",
            weights=[0.5, 0.25, 0.1, 0.15],
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

        #with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=4,
            allocation="manual",
            method="random",
            weights=[0.5, 0.25, 0.1, 0.15],
            mindist=90,
        ).samples_as_wkt())

        assert len(samples) > 90 #mindist means we may not get the full 100
        self.check_points_in_bounds(srast, samples)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

        #with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=4,
            allocation="manual",
            method="random",
            weights=[0.5, 0.25, 0.1, 0.15],
            access=self.access,
            layer_name='access',
            buff_inner=120,
            buff_outer=1200,
        ).samples_as_wkt())

        assert len(samples) == 100
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 120, 1200)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

        #with mindist and access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=4,
            allocation="manual",
            method="random",
            weights=[0.5, 0.25, 0.1, 0.15],
            access=self.access,
            layer_name='access',
            buff_outer=1200,
            mindist=90
        ).samples_as_wkt())

        assert len(samples) > 90 #mindist means we may not get the full 90 
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 0, 1200)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

    def test_random_allocation_optim(self):
        #test with mrast band = 0
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 10})
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=10,
            allocation="optim",
            mrast=self.rast,
            mrast_band=0,
            method="random",
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        
        #these percentages have already been pre-calculated
        assert percentages[0] - .18 == pytest.approx(0, abs=0.02)
        assert percentages[1] - .12 == pytest.approx(0, abs=0.02)
        assert percentages[2] - .10 == pytest.approx(0, abs=0.02)
        assert percentages[3] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[4] - .06 == pytest.approx(0, abs=0.02)
        assert percentages[5] - .05 == pytest.approx(0, abs=0.02)
        assert percentages[6] - .05 == pytest.approx(0, abs=0.02)
        assert percentages[7] - .04 == pytest.approx(0, abs=0.02)
        assert percentages[8] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[9] - .26 == pytest.approx(0, abs=0.02)
       
        #test with mrast band = 1
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=10,
            allocation="optim",
            mrast=self.rast,
            mrast_band=1,
            method="random",
        ).samples_as_wkt())
    
        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        
        #these percentages have already been pre-calculated
        assert percentages[0] - .12 == pytest.approx(0, abs=0.02)
        assert percentages[1] - .15 == pytest.approx(0, abs=0.02)
        assert percentages[2] - .14 == pytest.approx(0, abs=0.02)
        assert percentages[3] - .11 == pytest.approx(0, abs=0.02)
        assert percentages[4] - .10 == pytest.approx(0, abs=0.02)
        assert percentages[5] - .09 == pytest.approx(0, abs=0.02)
        assert percentages[6] - .08 == pytest.approx(0, abs=0.02)
        assert percentages[7] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[8] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[9] - .07 == pytest.approx(0, abs=0.02)

    def test_queinnec_allocation_equal(self):
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 5})

        #without mindist or access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast, 
            'strat_zq90',
            num_samples=500,
            num_strata=5,
            allocation="equal",
            method="Queinnec",
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0)

        #with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=5,
            allocation="equal",
            method="Queinnec",
            mindist=150,
        ).samples_as_wkt())
       
        assert len(samples) > 490 #mindist means we may not get the full 500 
        self.check_points_in_bounds(srast, samples)
        self.check_mindist(150, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0, abs=0.03)

        #with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=5,
            allocation="equal",
            method = "Queinnec",
            access = self.access,
            layer_name = 'access',
            buff_inner = 90,
            buff_outer = 300,
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 90, 300)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0, abs=0.03)

        #with access and mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=5,
            allocation="equal",
            method="Queinnec",
            access = self.access,
            layer_name = 'access',
            buff_outer=600,
            mindist=90,
        ).samples_as_wkt())

        assert len(samples) > 90 #mindist means we may not get the full 90 
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 0, 600)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.2 == pytest.approx(0, abs=0.03)

    def test_queinnec_allocation_proportional(self):
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 8})

        #without mindist or access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast, 
            'strat_zq90',
            num_samples=500,
            num_strata=8,
            allocation="prop",
            method="Queinnec"
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

        #with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=8,
            allocation="prop",
            method="Queinnec",
            mindist=150,
        ).samples_as_wkt())
        
        assert len(samples) > 490 #mindist means we may not get the full 500 
        self.check_points_in_bounds(srast, samples)
        self.check_mindist(150, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

        #with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=8,
            allocation="prop",
            method = "Queinnec",
            access = self.access,
            layer_name = 'access',
            buff_inner = 90,
            buff_outer = 300,
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 90, 300)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

        #with access and mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=128,
            num_strata=8,
            allocation="prop",
            method="Queinnec",
            access = self.access,
            layer_name = 'access',
            buff_outer=600,
            mindist=90,
        ).samples_as_wkt())

        assert len(samples) > 115 #mindist means we may not get the full 128 
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 0, 600)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        for percentage in percentages.values():
            assert percentage - 0.125 == pytest.approx(0, abs=0.03)

    def test_queinnec_allocation_manual(self):
        #without mindist or access
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 4})
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=4,
            allocation="manual",
            method="Queinnec",
            weights=[0.5, 0.25, 0.1, 0.15],
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

        #with mindist
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=4,
            allocation="manual",
            method="Queinnec",
            weights=[0.5, 0.25, 0.1, 0.15],
            mindist=90,
        ).samples_as_wkt())

        assert len(samples) > 90 #mindist means we may not get the full 100 
        self.check_points_in_bounds(srast, samples)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

        #with access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=4,
            allocation="manual",
            method="Queinnec",
            weights=[0.5, 0.25, 0.1, 0.15],
            access=self.access,
            layer_name='access',
            buff_inner=120,
            buff_outer=1200,
        ).samples_as_wkt())

        assert len(samples) == 100
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 120, 1200)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

        #with mindist and access
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=4,
            allocation="manual",
            method="Queinnec",
            weights=[0.5, 0.25, 0.1, 0.15],
            access=self.access,
            layer_name='access',
            buff_outer=1200,
            mindist=90
        ).samples_as_wkt())

        assert len(samples) > 90 #mindist means we may not get the full 100 
        self.check_points_in_bounds(srast, samples)
        self.check_access(samples, 0, 1200)
        self.check_mindist(90, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        assert percentages[0] - 0.5 == pytest.approx(0, abs=0.03)
        assert percentages[1] - 0.25 == pytest.approx(0, abs=0.03)
        assert percentages[2] - 0.1 == pytest.approx(0, abs=0.03)
        assert percentages[3] - 0.15 == pytest.approx(0, abs=0.03)

    def test_queinnec_allocation_optim(self):
        #test with mrast band = 0
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 10})
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=10,
            allocation="optim",
            mrast=self.rast,
            mrast_band=0,
            method="Queinnec",
        ).samples_as_wkt())

        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        
        #these percentages have already been pre-calculated
        assert percentages[0] - .18 == pytest.approx(0, abs=0.02)
        assert percentages[1] - .12 == pytest.approx(0, abs=0.02)
        assert percentages[2] - .10 == pytest.approx(0, abs=0.02)
        assert percentages[3] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[4] - .06 == pytest.approx(0, abs=0.02)
        assert percentages[5] - .05 == pytest.approx(0, abs=0.02)
        assert percentages[6] - .05 == pytest.approx(0, abs=0.02)
        assert percentages[7] - .04 == pytest.approx(0, abs=0.02)
        assert percentages[8] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[9] - .26 == pytest.approx(0, abs=0.02)
       
        #test with mrast band = 1
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=10,
            allocation="optim",
            mrast=self.rast,
            mrast_band=1,
            method="Queinnec",
        ).samples_as_wkt())
    
        assert len(samples) == 500
        self.check_points_in_bounds(srast, samples)
        percentages = self.get_allocation_percentages(srast, samples)
        
        #these percentages have already been pre-calculated
        assert percentages[0] - .12 == pytest.approx(0, abs=0.02)
        assert percentages[1] - .15 == pytest.approx(0, abs=0.02)
        assert percentages[2] - .14 == pytest.approx(0, abs=0.02)
        assert percentages[3] - .11 == pytest.approx(0, abs=0.02)
        assert percentages[4] - .10 == pytest.approx(0, abs=0.02)
        assert percentages[5] - .09 == pytest.approx(0, abs=0.02)
        assert percentages[6] - .08 == pytest.approx(0, abs=0.02)
        assert percentages[7] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[8] - .07 == pytest.approx(0, abs=0.02)
        assert percentages[9] - .07 == pytest.approx(0, abs=0.02)

    def test_existing_samples(self):
        #with force=True, test a couple different combinations of input parameters
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 10})
        existing = gpd.read_file(existing_shapefile_path)['geometry']
        existing_sample_count = len(existing)

        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=10,
            existing=self.existing,
            force=True,
            allocation="optim",
            mrast=self.rast,
            mrast_band=0,
            method="random"
        ).samples_as_wkt())

        for sample in existing:
            assert samples.contains(sample).any()

        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=10,
            existing=self.existing,
            force=True,
            allocation="optim",
            mrast=self.rast,
            mrast_band=0,
            method="Queinnec"
        ).samples_as_wkt())

        for sample in existing:
            assert samples.contains(sample).any()

        #with force=False, test a couple different combinations of input parameters
        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=100,
            num_strata=10,
            existing=self.existing,
            force=False,
            allocation="optim",
            mrast=self.rast,
            mrast_band=0,
            method="random"
        ).samples_as_wkt())
    
        existing_in_final = 0
        for sample in existing:
            if samples.contains(sample).any():
                existing_in_final += 1

        assert existing_in_final != 0
        assert existing_in_final != existing_sample_count

        samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
            srast,
            'strat_zq90',
            num_samples=500,
            num_strata=10,
            existing=self.existing,
            force=False,
            allocation="manual",
            weights=[0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, .01],
            mrast=self.rast,
            mrast_band=0,
            method="Queinnec"
        ).samples_as_wkt())

        existing_in_final = 0
        for sample in existing:
            if samples.contains(sample).any():
                existing_in_final += 1

        assert existing_in_final != 0
        assert existing_in_final != existing_sample_count

    def test_queinnec_focal_window(self):
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 5})

        #queinnec sampling works by first adding pixels which have a focal
        #window containing the same pixels, then automatically moving to 
        #random sampling if required. To test the focal window calculation
        #adequately while not forcing the algorithm to move to random sampling,
        #sampling a small amount of pixels, but doing it repeatedly.

        #test 3x3 focal window
        for _ in range(20):
            samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
                srast,
                'strat_zq90',
                num_samples=50,
                num_strata=5,
                wrow=3,
                wcol=3,
                allocation="equal",
                method="Queinnec",
            ).samples_as_wkt())
            self.check_focal_window(srast, samples, 3, 3)

        #test 5x5 focal window
        for _ in range(100):
            samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
                srast,
                'strat_zq90',
                num_samples=5,
                num_strata=5,
                wrow=5,
                wcol=5,
                allocation="equal",
                method="Queinnec",
            ).samples_as_wkt())
            self.check_focal_window(srast, samples, 5, 5)

        #test 3x5 focal window
        for _ in range(100):
            samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
                srast,
                'strat_zq90',
                num_samples=5,
                num_strata=5,
                wrow=3,
                wcol=5,
                allocation="equal",
                method="Queinnec",
            ).samples_as_wkt())
            self.check_focal_window(srast, samples, wrow=3, wcol=5)

        #test 5x3 focal window
        for _ in range(100):
            samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(
                srast,
                'strat_zq90',
                num_samples=5,
                num_strata=5,
                wrow=5,
                wcol=3,
                allocation="equal",
                method="Queinnec",
            ).samples_as_wkt())
            self.check_focal_window(srast, samples, wrow=5, wcol=3)
 
    def test_function_inputs(self):
        srast = sgs.stratify.quantiles(self.rast, num_strata={"zq90": 5})

        #test wrow inputs
        for wrow in [-1, 0, 2]:
            with pytest.raises(ValueError):
                sgs.sample.strat(
                    srast,
                    'strat_zq90',
                    num_strata=5,
                    num_samples=5,
                    wrow=wrow,
                    wcol=3,
                    allocation="equal",
                    method="Queinnec",
                )

        #test wcol inputs
        for wcol in [-1, 0, 2]:
            with pytest.raises(ValueError):
                sgs.sample.strat(
                    srast,
                    'strat_zq90',
                    num_strata=5,
                    num_samples=5,
                    wrow=3,
                    wcol=wcol,
                    allocation="equal",
                    method="Queinnec",
                )


        #test method strings
        for method in ["queinnec", "Random", "", "badstring"]:
            with pytest.raises(ValueError):
                sgs.sample.strat(
                    srast,
                    'strat_zq90',
                    num_strata=5,
                    num_samples=5,
                    allocation="equal",
                    method=method,
                )

        #test allocation strings
        for allocation in ["", "badstring", "Equal", "propp"]:
            with pytest.raises(ValueError):
                sgs.sample.strat(
                    srast,
                    'strat_zq90',
                    num_strata=5,
                    num_samples=5,
                    allocation=allocation,
                    method="random",
                )

        #test weights
        for weights in [None, [], [0.2, 0.2, 0.2, 0.2, 0.21], [0.5, 0.4, 0.05, 0.05]]:
            with pytest.raises(ValueError):
                sgs.sample.strat(
                    srast,
                    'strat_zq90',
                    num_strata=5,
                    num_samples=5,
                    allocation="manual",
                    weights=weights,
                    method="random",
                )

        #test mindist
        with pytest.raises(ValueError):
            sgs.sample.strat(
                srast,
                'strat_zq90',
                num_strata=5,
                num_samples=5,
                allocation="equal",
                method="random",
                mindist=-1,
            )
