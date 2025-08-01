import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster_small.tif')
samples = sgs.sample.balanced(rast, algorithm="lcube", num_samples=200, plot=True)

