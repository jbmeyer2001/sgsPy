import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
print('lcube')
samples = sgs.sample.balanced(rast, algorithm="lcube", num_samples=200, plot=True)
print('lpm2_kdtree')
samples = sgs.sample.balanced(rast, algorithm='lpm2_kdtree', num_samples=200, plot=True)
print('lcubestratified')
srast = sgs.stratify.quantiles(rast, num_strata={'pzabove2':5})
samples = sgs.sample.balanced(rast, algorithm='lcubestratified', srast=srast, srast_band='strat_pzabove2', num_samples=200, plot=True)
