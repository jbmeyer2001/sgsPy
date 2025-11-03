import sgs
import sys
#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/RMF_DEM1m.tif')
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
existing = sgs.SpatialVector('/home/jbmeyer/extdata/existing.shp')
#access = sgs.SpatialVector('/home/jbmeyer/RMF_vector_files/RMF_Roads_clip.shp')
access = sgs.SpatialVector('/home/jbmeyer/extdata/access.shp')

#print("creating DEM_quantiles.tif")
#srast = sgs.stratify.quantiles(rast, num_strata=5, filename="out/DEM_quantiles.tif")
#print("DEM_quantiles.tif created")
#srast = sgs.SpatialRaster("out/DEM_quantiles.tif")
srast = sgs.stratify.quantiles(rast, num_strata={'pzabove2': 5})

#samples =sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, allocation="equal", method="Queinnec", filename="out/samples.shp")
samples =sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, allocation="equal", method="Queinnec", plot=True)
print("finished sampling no access no existing")

#samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, access=access, layer_name="RMF_Roads_clip", buff_outer=1000, allocation="equal", method="Queinnec", filename="out/samples_access.shp")
samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, access=access, layer_name="access", buff_outer=120, allocation="equal", method="Queinnec", plot=True)
print("finished sampling access")

#samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata=5, existing=existing, allocation="equal", method="Queinnec", filename="out/samples_existing.shp")
samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata=5, existing=existing, allocation="equal", method="Queinnec", plot=True)
print("finished sampling existing")

