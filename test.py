import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster_small.tif')
print('lcube')
samples = sgs.sample.balanced(rast, algorithm="lcube", num_samples=50, plot=True)
print()
print('lpm2_kdtree')
samples = sgs.sample.balanced(rast, algorithm='lpm2_kdtree', num_samples=50, plot=True)
print()
print('lcubestratified')
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
srast = sgs.stratify.quantiles(rast, num_strata={'pzabove2':5})
samples = sgs.sample.balanced(rast, algorithm='lcubestratified', srast=srast, srast_band='strat_pzabove2', num_samples=200, plot=True)
