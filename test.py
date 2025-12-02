import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
print("HERE 1")
samples = sgs.sample.clhs(rast, num_samples=200, plot=True)
print("HERE 3")

