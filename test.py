import sgs
import matplotlib.pyplot as plt

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')

samples = sgs.sample.systematic(
    rast,
    500,
    'square',
    'centers',
    False
)
