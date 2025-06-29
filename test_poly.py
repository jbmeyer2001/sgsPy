import sgs

vect = sgs.SpatialVector('/home/jbmeyer/extdata/inventory_polygons.shp')
rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')

strat_rast = sgs.poly(rast, vect, 'NUTRIENTS', 'inventory_polygons', ['rich', 'medium', 'poor']) 
strat_rast.plot()
