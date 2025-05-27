import sgs
from sgs.utils import SpatialRaster

rast = SpatialRaster('/home/jbmeyer/extdata/mraster_small.tif')
sgs.balanced(rast, 5)
