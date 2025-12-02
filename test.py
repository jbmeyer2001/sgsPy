import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
print("HERE 1")
samples = sgs.sample.clhs(rast, num_samples=10, plot=True)
print("HERE 3")

