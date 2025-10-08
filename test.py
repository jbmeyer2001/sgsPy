import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
pcomp = sgs.calculate.pca(rast, 5)
