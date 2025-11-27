import sgs

print("HERE 1")
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
print("HERE 2")
srast = sgs.stratify.quantiles(rast, num_strata={rast.bands[0]:10}, filename="test.tif")
print("HERE 3")
srast.plot(band=0)

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Merged_Aspect_DEM5m.tif')
#pcomp = sgs.calculate.pca(rast, 1, filename='aspect_dem_pca.tif')

