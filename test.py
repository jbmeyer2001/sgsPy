import sgs
import rasterio
import matplotlib.pyplot as plt

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
