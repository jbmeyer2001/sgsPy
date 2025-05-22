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

__all__ = [
    'access_shapefile_path',
    'existing_shapefile_path',
    'existing_geodatabase_path',
    'existing_geojson_path',
    'inventory_polygons_shapefile_path',
    'mraster_geotiff_path',
    'mraster_small_geotiff_path',
    'sraster_geotiff_path',
]
