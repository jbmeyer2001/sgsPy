import sgs

print("HERE 1")
rast = sgs.SpatialRaster('E:/test_delete_later/merged.tif')
print("HERE 2")
pcomp = sgs.calculate.pca(rast, 3, filename="E:/test_delete_later/large_pca_out.tif")
print("HERE 3")
pcomp.plot(band=0)
pcomp.plot(band=1)
pcomp.plot(band=2)

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Merged_Aspect_DEM5m.tif')
#pcomp = sgs.calculate.pca(rast, 1, filename='aspect_dem_pca.tif')

