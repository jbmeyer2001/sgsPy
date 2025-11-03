import sgs
import sys
rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/RMF_DEM1m.tif')
existing = sgs.SpatialVector('/home/jbmeyer/extdata/existing.shp')
access = sgs.SpatialVector('/home/jbmeyer/RMF_vector_files/RMF_Roads_clip.shp')

#print("creating DEM_quantiles.tif")
#srast = sgs.stratify.quantiles(rast, num_strata=5, filename="out/DEM_quantiles.tif")
#print("DEM_quantiles.tif created")
srast = sgs.SpatialRaster("out/DEM_quantiles.tif")

samples =sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, allocation="equal", method="random", filename="out/samples.shp")
print("finished sampling no access no existing")

samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, access=access, layer_name="RMF_Roads_clip", buff_outer=1000, allocation="equal", method="random", filename="out/samples_access.shp")
print("finished sampling access")

samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata=5, existing=existing, allocation="equal", method="random", filename="out/samples_existing.shp")
print("finished sampling existing")

