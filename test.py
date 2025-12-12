import sgs
import rasterio
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
srast = sgs.stratify.quantiles(rast, num_strata={"zq90":10, "pzabove2":10})
srast_rasterio = srast.to_rasterio()
