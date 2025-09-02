
import time
import sgs
#rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
#srast = sgs.stratify.breaks(rast, breaks={"zq90": [3, 5, 7, 9, 11, 13, 15]})
#srast.plot()

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/RMF_DEM1m.tif')
#start = time.time()
#srast = sgs.stratify.breaks(rast, list(range(10, 490, 10)), filename="test_large_image_breaks.tif")
#end = time.time()
#print(str(end - start))

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/CHM.tif')
#start = time.time()
#srast = sgs.stratify.breaks(rast, list(range(1, 49, 1)), filename="strat_CHM.tif")
#end = time.time()
#print(str(end - start))

rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/CHM_DEM.tif')
start = time.time()
srast = sgs.stratify.breaks(rast, [list(range(10, 490, 10)), list(range(1, 49, 1))], map=True)
end = time.time()
print(str(end - start))
