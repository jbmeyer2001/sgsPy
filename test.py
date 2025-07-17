import sgs
import matplotlib.pyplot as plt

rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')

for location in ['centers', 'corners', 'random']:
    for shape in ['square', 'hexagon']:
        samples = sgs.sample.systematic(
            rast,
            500,
            shape,
            location,
            plot=True
        )
        #print(samples)

