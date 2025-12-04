import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Merged_Aspect_DEM5m.tif')
print("HERE!")
samples = sgs.sample.clhs(rast, num_samples=200, plot=True)

