import sgs
import sys

import geopandas as gpd

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/CHM.tif')
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
existing = sgs.SpatialVector('/home/jbmeyer/extdata/existing.shp')
#access = sgs.SpatialVector('/home/jbmeyer/RMF_vector_files/RMF_Roads_clip.shp')
access = sgs.SpatialVector('/home/jbmeyer/extdata/access.shp')

#print("creating CHM_quantiles.tif")
#srast = sgs.stratify.quantiles(rast, num_strata=10, filename="out/CHM_quantiles.tif")
#srast = sgs.SpatialRaster("out/CHM_quantiles.tif")
#print("CHM_quantiles.tif created")
#srast = sgs.SpatialRaster("out/DEM_quantiles.tif")
srast = sgs.stratify.quantiles(rast, num_strata={'zq90': 10})

#samples =sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, allocation="equal", method="Queinnec", filename="out/samples.shp")
#samples =sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, allocation="equal", method="Queinnec", plot=True)
#print("finished sampling no access no existing")

#samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, access=access, layer_name="RMF_Roads_clip", buff_outer=1000, allocation="equal", method="Queinnec", filename="out/samples_access.shp")
#samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata = 5, access=access, layer_name="access", buff_outer=120, allocation="equal", method="Queinnec", plot=True)
#print("finished sampling access")

#samples = sgs.sample.strat(srast, 0, num_samples=100, num_strata=10, existing=existing, allocation="equal", method="Queinnec", filename="out/samples_existing.shp")

#samples = sgs.sample.strat(srast, 0, num_samples=100, num_strata=10, existing=existing, force=True, allocation="equal", method="Queinnec", filename="out/samples_existing_forced.shp")

#samples = sgs.sample.strat(srast, 0, num_samples=500, num_strata=5, existing=existing, allocation="equal", method="Queinnec", plot=True)
#print("finished sampling existing")

#samples = sgs.sample.strat(srast, 0, num_samples = 500, num_strata=10, mrast=rast, mrast_band=0, allocation="optim", method="Queinnec", filename="out/samples_queinnec_optim.shp")





#existing_samples = gpd.read_file('/home/jbmeyer/extdata/existing.shp')['geometry']
samples = gpd.GeoSeries.from_wkt(sgs.sample.strat(srast, 'strat_zq90', num_samples = 100, num_strata = 10, existing=existing, force=True, allocation="optim", mrast=rast, mrast_band=0, method="random", plot=True).samples_as_wkt())

print()
print()
print()
print()
print()

samples = gpd.GeoSeries.from_wkt(sgs.sample.srs(rast, num_samples=100, existing=existing, plot=True).samples_as_wkt())

#for sample in samples:
#    print (sample)

print()
print()
print()
print()
print()

for sample in existing_samples:
    plot(sample)

#for sample in existing_samples:
#    if not samples.contains(sample).any():
#        print(sample)



