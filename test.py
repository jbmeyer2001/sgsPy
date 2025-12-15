import sgs
import geopandas as gpd
import matplotlib.pyplot as plt

"""
rast = sgs.SpatialRaster("C:/Users/jmeyer03/projects/Github/sgs/tests/files/mraster.tif")
srast = sgs.stratify.quantiles(rast, num_strata={"zq90":5})
ds, arr = srast.to_gdal(with_arr=True)
arr[arr == 2] = 4
print('arr[0]')
print(arr[0])
print()
srast = sgs.SpatialRaster.from_gdal(ds, arr)
print('srast.band(0)')
print(srast.band(0))
sample = sgs.sample.strat(strat_rast=srast, band=0, num_samples=200, num_strata=5, plot=True)

rast = sgs.SpatialRaster("C:/Users/jmeyer03/projects/Github/sgs/tests/files/mraster.tif")
srast = sgs.stratify.quantiles(rast, num_strata={"zq90":5, "pzabove2":5})
srasterio, arr = srast.to_rasterio(with_arr = True)
arr[arr == 2] = 4
print('arr[0]')
print(arr[0])
print()
print('arr[1]')
print(arr[1])
print()
srast = sgs.SpatialRaster.from_rasterio(srasterio, arr)
print('srast.band(0)')
print(srast.band(0))
print()
print("srast.band(1)")
print(srast.band(1))
print()
sample = sgs.sample.strat(strat_rast = srast, band=0, num_samples=200, num_strata = 5, plot=True)
"""


rast = sgs.SpatialRaster("/home/jbmeyer/extdata/mraster.tif")
#gdf = gpd.read_file("/home/jbmeyer/extdata/access.shp")
#access = sgs.SpatialVector.from_geopandas(gdf)
#samples = sgs.sample.srs(rast=rast, num_samples=250, access=access, buff_outer=500 , plot=True) 

gdf = gpd.read_file("/home/jbmeyer/extdata/inventory_polygons.shp")
print("HURR 1")
polygons = sgs.SpatialVector.from_geopandas(gdf)
print("HURR 2")
srast = sgs.stratify.poly(rast, polygons, attribute='NUTRIENTS', layer_name='inventory_polygons', features=['poor', 'medium', 'rich'])
print("HURR 3")
srast.plot()
