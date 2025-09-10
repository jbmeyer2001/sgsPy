# python ./test_large_CHM.py ../../RMF_LiDAR_Metrics/CHM.tif CHM_test.tif

import sys
import numpy as np
import rasterio
from rasterio.windows import Window

image_width = 106000
image_height = 118000

rast_name = sys.argv[1]
srast_name = sys.argv[2]

with rasterio.open(sys.argv[1]) as rast:
    with rasterio.open(sys.argv[2]) as srast:
        for y in range(image_height):
            rast_scanline = rast.read(1, window=Window(0, y, image_width, 1))
            srast_scanline = srast.read(1, window=Window(0, y, image_width, 1))

            for x in range(image_width):
                val = rast_scanline[0, x]
                strata = srast_scanline[0, x]

                if val < 0:
                    if strata > -1:
                        msg = "[NAN PROBLEM] value is " + str(rast_scanline[x]) + " but strata is " + str(strata) + "at index: [" + str(x) + "][" + str(y) + "]" 
                        raise RuntimeError(msg)
                    continue

                val = np.trunc(val - .00001)
                if val < 0 : val = 0
                if val > 48 : val = 48

                if val != strata:
                    msg = "value is " + str(rast_scanline[0, x]) + " but strata is " + str(strata) + " at index: [" + str(x) + "][" + str(y) + "]" 
                    raise RuntimeError(msg)

            if y % 100 == 0:
                print(y)
