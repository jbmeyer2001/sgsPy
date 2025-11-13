import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/merged.tif')
pcomp = sgs.calculate.pca(rast, 3, filename="large_pca_out.tif")
#pcomp.plot(band=0)
#pcomp.plot(band=1)
#pcomp.plot(band=2)

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Merged_Aspect_DEM5m.tif')
#pcomp = sgs.calculate.pca(rast, 1, filename='aspect_dem_pca.tif')

