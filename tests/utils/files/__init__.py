import os

folder = os.path.dirname(__file__)

access_shapefile_path = os.path.join(folder, 'access.shp')
existing_shapefile_path = os.path.join(folder, 'existing.shp')
existing_geodatabase_path = os.path.join(folder, 'existing.gdb')
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

__all__ = [
    'access_shapefile_path',
    'existing_shapefile_path',
    'existing_geodatabase_path',
    'existing_geojson_path',
    'inventory_polygons_shapefile_path',
    'mraster_geotiff_path',
    'mraster_small_geotiff_path',
    'sraster_geotiff_path',
    'sraster2_geotiff_path',
    'mraster_zq90_path',
    'mraster_pzabove2_path',
    'mraster_zsd_path',
    'mraster_small_zq90_path',
    'mraster_small_pzabove2_path',
    'mraster_small_zsd_path',
    'sraster_strata_path',
    'sraster2_band_path',
]
