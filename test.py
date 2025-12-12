import sgs
from osgeo import gdal
import matplotlib.pyplot as plt

ds = gdal.Open("C:/Users/jmeyer03/projects/Github/sgs/tests/files/mraster.tif")
rast = sgs.SpatialRaster.from_gdal(ds)
srast = sgs.stratify.quantiles(rast, num_strata={rast.bands[0]: 10})
srast_gdal = srast.to_gdal()

arr = srast_gdal.GetRasterBand(1).ReadAsArray()
plt.imshow(arr)
plt.show()

"""
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
srast = sgs.stratify.quantiles(rast, num_strata={"zq90":5})
srasterio, arr = srast.to_rasterio(with_arr=True)
arr[arr == 2] = 4
print('arr[0]')
print(arr[0])
print()
srast = sgs.SpatialRaster.from_rasterio(srasterio, arr)
print('srast.band(0)')
print(srast.band(0))
sample = sgs.sample.strat(strat_rast=srast, band=0, num_samples=200, num_strata=5, plot=True)
"""
