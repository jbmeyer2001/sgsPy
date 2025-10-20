import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
pcomp = sgs.calculate.pca(rast, 3)
pcomp.plot(band=0)
pcomp.plot(band=1)
pcomp.plot(band=2)
