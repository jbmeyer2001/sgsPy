import sgs
from sgs.utils import SpatialRaster

print("imported SpatialRaster!")

rast = SpatialRaster('/home/jbmeyer/extdata/mraster.tif')

print("opened mraster.tif")

rast.load_arr()

print("loaded array")
