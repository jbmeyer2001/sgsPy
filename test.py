import sgs
from sgs.utils import SpatialRaster

rast = SpatialRaster('/home/jbmeyer/extdata/sraster.tif')
rast.info()
rast.plot_image()
