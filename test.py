import sgs

rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/RMF_DEM1m.tif')
#rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
access = sgs.SpatialVector('/home/jbmeyer/RMF_vector_files/RMF_Roads_clip.shp')
#access = sgs.SpatialVector('/home/jbmeyer/extdata/access.shp')
existing = sgs.SpatialVector('/home/jbmeyer/extdata/existing.shp')

print("square centers")
samples = sgs.sample.systematic(rast, 10000, "square", "centers", filename='out/square_centers.shp')
print("square corners")
samples = sgs.sample.systematic(rast, 10000, "square", "corners", filename='out/square_corners.shp')
print("square random")
samples = sgs.sample.systematic(rast, 10000, "square", "random", filename='out/square_random.shp')
print("hexagon centers")
samples = sgs.sample.systematic(rast, 10000, "hexagon", "centers", filename='out/hexagon_centers.shp')
print("hexagon corners")
samples = sgs.sample.systematic(rast, 10000, "hexagon", "corners", filename='out/hexagon_corners.shp')
print("hexagon random")
samples = sgs.sample.systematic(rast, 10000, "hexagon", "random", filename='out/hexagon_random.shp')

"""
Testing w/ force

samples = sgs.sample.systematic(rast, 10000, "square", "centers", force=True, filename='out/square_centers.shp')
samples = sgs.sample.systematic(rast, 10000, "square", "corners", force=True, filename='out/square_corners.shp')
samples = sgs.sample.systematic(rast, 10000, "square", "random", force=True, filename='out/square_random.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "centers", force=True, filename='out/hexagon_centers.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "corners", force=True, filename='out/hexagon_corners.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "random", force=True, filename='out/hexagon_random.shp')
"""

"""
Testing w/ access

samples = sgs.sample.systematic(rast, 10000, "square", "centers", access=access, layer_name='RMF_Roads_clip', buff_outer=500, filename='out/square_centers.shp')
samples = sgs.sample.systematic(rast, 10000, "square", "corners", access=access, layer_name='RMF_Roads_clip', buff_outer=500, filename='out/square_corners.shp')
samples = sgs.sample.systematic(rast, 10000, "square", "random", access=access, layer_name='RMF_Roads_clip', buff_outer=500, filename='out/square_random.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "centers", access=access, layer_name='RMF_Roads_clip', buff_outer=500, filename='out/hexagon_centers.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "corners", access=access, layer_name='RMF_Roads_clip', buff_outer=500, filename='out/hexagon_corners.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "random", access=access, layer_name='RMF_Roads_clip', buff_outer=500, filename='out/hexagon_random.shp')
"""

"""
Testing w/ existing

samples = sgs.sample.systematic(rast, 10000, "square", "centers", existing=existing, filename='out/square_centers.shp')
samples = sgs.sample.systematic(rast, 10000, "square", "corners", existing=existing, filename='out/square_corners.shp')
samples = sgs.sample.systematic(rast, 10000, "square", "random", existing=existing, filename='out/square_random.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "centers", existing=existing, filename='out/hexagon_centers.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "corners", existing=existing, filename='out/hexagon_corners.shp')
samples = sgs.sample.systematic(rast, 10000, "hexagon", "random", existing=existing, filename='out/hexagon_random.shp')
"""


