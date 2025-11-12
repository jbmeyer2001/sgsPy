import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
pcomp = sgs.calculate.pca(rast, 3, filename='sgspy_pca_result.tif')
pcomp.plot(band=1)
#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Merged_Aspect_DEM5m.tif')
#pcomp = sgs.calculate.pca(rast, 1, filename='aspect_dem_pca.tif')

