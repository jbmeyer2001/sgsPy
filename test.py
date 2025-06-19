import os

import sgs
from sgs import SpatialRaster
import numpy as np

folder = os.path.join(os.path.dirname(__file__), 'tests', 'utils', 'files')
access_shapefile_path = os.path.join(folder, 'access.shp')
existing_shapefile_path = os.path.join(folder, 'existing.shp')
existing_geojson_path = os.path.join(folder, 'existing.geojson')
inventory_polygons_shapefile_path = os.path.join(folder, 'inventory_polygons.shp')
mraster_geotiff_path = os.path.join(folder, 'mraster.tif')
mraster_small_geotiff_path = os.path.join(folder, 'mraster_small.tif')
sraster_geotiff_path = os.path.join(folder, 'sraster.tif')
sraster2_geotiff_path = os.path.join(folder, 'sraster2.tif')
mraster_zq90_path = os.path.join(folder, 'mraster_zq90.npy')
mraster_pzabove2_path = os.path.join(folder, 'mraster_pzabove2.npy')
mraster_zsd_path = os.path.join(folder, 'mraster_zsd.npy')
mraster_small_zq90_path = os.path.join(folder, 'mraster_small_zq90.npy')
mraster_small_pzabove2_path = os.path.join(folder, 'mraster_small_pzabove2.npy')
mraster_small_zsd_path = os.path.join(folder, 'mraster_small_zsd.npy')
sraster_strata_path = os.path.join(folder, 'sraster_strata.npy')
sraster2_band_path = os.path.join(folder, 'sraster2_band.npy')

rast = SpatialRaster(mraster_geotiff_path)
assert rast.width == 373
assert rast.height == 277
assert rast.band_count == 3
assert rast.pixel_width == 20
assert rast.pixel_height == 20
assert rast.xmin == 431100
assert rast.xmax == 438560
assert rast.ymin == 5337700
assert rast.ymax == 5343240
assert np.array_equal(np.load(mraster_zq90_path), rast[0], equal_nan=True)
assert np.array_equal(np.load(mraster_pzabove2_path), rast[1], equal_nan=True)
assert np.array_equal(np.load(mraster_zsd_path), rast[2], equal_nan=True)

rast = SpatialRaster(mraster_small_geotiff_path)#    def mraster_small_check(self, rast):
assert rast.width == 141
assert rast.height == 110
assert rast.band_count == 3
assert rast.pixel_width == 20
assert rast.pixel_height == 20
assert rast.xmin == 432440
assert rast.xmax == 435260
assert rast.ymin == 5340000
assert rast.ymax == 5342200
assert np.array_equal(np.load(mraster_small_zq90_path), rast[0], equal_nan=True)
assert np.array_equal(np.load(mraster_small_pzabove2_path), rast[1], equal_nan=True)
assert np.array_equal(np.load(mraster_small_zsd_path), rast[2], equal_nan=True)
    
rast = SpatialRaster(sraster_geotiff_path)#def sraster_check(self, rast):
assert rast.width == 373
assert rast.height == 277
assert rast.band_count == 1
assert rast.pixel_width == 20
assert rast.pixel_height == 20
assert rast.xmin == 431100
assert rast.xmax == 438560
assert rast.ymin == 5337700
assert rast.ymax == 5343240
assert np.array_equal(np.load(sraster_strata_path), rast[0], equal_nan=True)

rast = SpatialRaster(sraster2_geotiff_path)#    def sraster2_check(self, rast):
assert rast.width == 861
assert rast.height == 611
assert rast.band_count == 1
assert abs(rast.pixel_width - 0.0003353943) < 1e-10# pytest.approx(0, abs=1e-10)
assert abs(rast.pixel_height - 0.0003353943) < 1e-10#== pytest.approx(0, abs=1e-10)
assert abs(rast.xmin - -71.9365) < 1e-4#== pytest.approx(0, abs=1e-4)
assert abs(rast.xmax - -71.64773) < 1e-4#== pytest.approx(0, abs=1e-4)
assert abs(rast.ymin - 45.4883) < 1e-4#== pytest.approx(0, abs=1e-4)
assert abs(rast.ymax - 45.69332) <1e-4#== pytest.approx(0, abs=1e-4)
assert np.array_equal(np.load(sraster2_band_path), rast[0], equal_nan=True)

rast = sgs.utils.raster.SpatialRaster(mraster_small_geotiff_path)
zq90 = np.load(mraster_small_zq90_path)
pzabove2 = np.load(mraster_small_pzabove2_path)
zsd = np.load(mraster_small_zsd_path)

assert np.array_equal(zq90, rast[0], equal_nan=True)
assert np.array_equal(pzabove2, rast[1], equal_nan=True)
assert np.array_equal(zsd, rast[2], equal_nan=True)

assert np.array_equal(zq90, rast['zq90'], equal_nan=True)
assert np.array_equal(pzabove2, rast['pzabove2'], equal_nan=True)
assert np.array_equal(zsd, rast['zsd'], equal_nan=True)

assert np.array_equal(zq90[0:10, 0:10], rast['zq90', 0:10, 0:10], equal_nan=True)
assert np.array_equal(zq90[131:141, 100:110], rast['zq90', 131:141, 100:110], equal_nan=True)
assert np.array_equal([zq90, pzabove2], rast[0:2], equal_nan=True)
assert np.array_equal(zq90[73, 46], rast['zq90', 73, 46], equal_nan=True)

vec = sgs.utils.vector.SpatialVector(access_shapefile_path)
assert len(vec.layers) == 1
assert vec.layers[0] == 'access'
layer_info = vec.cpp_vector.get_layer_info('access')
assert int(layer_info["feature_count"]) == 167
assert int(layer_info["field_count"]) == 2
assert layer_info["geometry_type"] == "Line String"
assert float(layer_info["xmin"]) == 431100
assert float(layer_info["xmax"]) == 438560
assert float(layer_info["ymin"]) == 5337700
assert float(layer_info["ymax"]) == 5343240

vec = sgs.utils.vector.SpatialVector(existing_shapefile_path)
assert len(vec.layers) == 1
assert vec.layers[0] == 'existing'
layer_info = vec.cpp_vector.get_layer_info('existing')
assert int(layer_info["feature_count"]) == 200
assert int(layer_info["field_count"]) == 1
assert layer_info["geometry_type"] == "Point"
assert float(layer_info["xmin"]) == 431110
assert float(layer_info["xmax"]) == 438530
assert float(layer_info["ymin"]) == 5337710
assert float(layer_info["ymax"]) == 5343230

vec = sgs.utils.vector.SpatialVector(inventory_polygons_shapefile_path)
assert len(vec.layers) == 1
assert vec.layers[0] == 'inventory_polygons'
layer_info = vec.cpp_vector.get_layer_info('inventory_polygons')
assert int(layer_info["feature_count"]) == 632
assert int(layer_info["field_count"]) == 3
assert layer_info["geometry_type"] == "Polygon"
assert float(layer_info["xmin"]) == 431100
assert float(layer_info["xmax"]) == 438560
assert float(layer_info["ymin"]) == 5337700
assert float(layer_info["ymax"]) == 5343240

vec = sgs.utils.vector.SpatialVector(existing_geojson_path)
assert len(vec.layers) == 1
assert vec.layers[0] == 'existing'
layer_info = vec.cpp_vector.get_layer_info('existing')
assert int(layer_info["feature_count"]) == 200
assert int(layer_info["field_count"]) == 1
assert layer_info["geometry_type"] == "Point"
assert float(layer_info["xmin"]) == 431110
assert float(layer_info["xmax"]) == 438530
assert float(layer_info["ymin"]) == 5337710
assert float(layer_info["ymax"]) == 5343230

print("PASSED raster.cpp and vector.cpp test")

#rast = SpatialRaster(sraster_geotiff_path)
#samples = sgs.srs(rast, mindist=200, num_samples=40, plot=True, filename="test_outputs/test_file_out.shp")
#print(samples)

#rast = SpatialRaster(mraster_small_geotiff_path)
#samples = sgs.srs(rast, num_samples=60, plot=True, filename="test_outputs/test_file_out.geojson")
#print(samples)

#rast = SpatialRaster(mraster_geotiff_path)
#samples = sgs.srs(rast, mindist=1000, num_samples=50, plot=True, filename="test_outputs/test_file_out.shp")
#print(samples)

rast = SpatialRaster(mraster_small_geotiff_path)
strat_rast = sgs.stratify.breaks(rast,[[10, 12], [50,80], [3,4]], map=True)
print(strat_rast[0])
print(strat_rast[1])
print(strat_rast[2])
print(strat_rast[3])
