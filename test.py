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

samples = sgs.sample.systematic(
    rast,
    500,
    'hexagon',
    'centers',
    False
)

samples = sgs.sample.systematic(
    rast,
    500,
    'triangle',
    'centers',
    False
)
