import sgs
import sys
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
#rast = sgs.SpatialRaster('home/jbmeyer/RMF_LiDAR_Metrics/CHM.tif')
srast = sgs.stratify.quantiles(rast, num_strata={"zq90": 5})
#srast = sgs.stratify.quantiles(rast, num_strata=5)

print(srast.bands)
print("about to start sampling using strat")
samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, allocation="equal", method="random", plot=True)
print("finished sampling!")

