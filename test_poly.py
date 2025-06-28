import sgs

vect = sgs.SpatialVector('/home/jbmeyer/extdata/inventory_polygons.shp')
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')

samples = sgs.poly(rast, vect, 'NUTRIENTS', ['rich', 'medium', 'poor']) 
