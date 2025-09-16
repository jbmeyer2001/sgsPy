import time
import sgs
import sys
#rast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
#srast = sgs.stratify.quantiles(rast, num_strata={"zq90": 5})
#srast.plot()
#srast = sgs.stratify.map((srast, ["strat_zq90", "strat_pzabove2"], [8, 4]), filename="test_map.tif")
#srast.plot()
#srast.plot(band=0)
#srast.plot(band=1)
#srast.plot(band=2)

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/RMF_DEM1m.tif')
#start = time.time()
#srast_dem = sgs.stratify.quantiles(rast, 10, eps=0.001, filename="srast_DEM.tif")
#end = time.time()
#print(str(end - start))

#rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/CHM.tif')
#start = time.time()
#srast_chm = sgs.stratify.quantiles(rast, 10, filename="srast_CHM")
#end = time.time()
#print(str(end - start))

rast = sgs.SpatialRaster('/home/jbmeyer/RMF_LiDAR_Metrics/Merged_Aspect_DEM5m.tif')
start = time.time()
srast = sgs.stratify.quantiles(rast, {rast.bands[1]: 10}, filename="5m_band1.tif", thread_count = 1)
#srast = sgs.stratify.breaks(rast, [list(range(10, 490, 10)), list(range(1, 49, 1))], map=True)
end = time.time()
print(str(end - start))
#srast.plot(band=0)
#srast.plot(band=1)
#srast.plot(band=2)


#srast_chm = sgs.SpatialRaster("srast_CHM.tif")
#srast_dem = sgs.SpatialRaster("srast_DEM.tif")
#start = time.time()
#srast = sgs.stratify.map((srast_chm, 0, 49), (srast_dem, 0, 49), filename="test_large_mapped.tif")
#end = time.time()
#print(str(end - start))
