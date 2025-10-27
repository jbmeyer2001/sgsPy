import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
access = sgs.SpatialVector('/home/jbmeyer/extdata/access.shp')
existing = sgs.SpatialVector('/home/jbmeyer/extdata/existing.shp')

samples = sgs.sample.srs(rast, num_samples = 150, plot=True)

