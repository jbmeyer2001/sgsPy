import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Zentropy.tif')
access = sgs.SpatialVector('/home/jbmeyer/RMF_vector_files/RMF_Roads_clip.shp')
existing = sgs.SpatialVector('/home/jbmeyer/extdata/existing.shp')

samples = sgs.sample.srs(rast, num_samples = 250, filename='out/samples.shp', existing=existing, access=access, layer_name='RMF_Roads_clip', buff_outer = 100)

#samples = sgs.sample.srs(rast, num_samples = 250, existing=existing, access=access, layer_name='access', buff_inner=30, buff_outer=180, mindist = 60, plot=True)

