# python ./test_large_CHM.py ../../RMF_LiDAR_Metrics/CHM.tif CHM_test.tif

import sys
import numpy as np
import rasterio
from rasterio.windows import Window

image_width = 106000
image_height = 118000

with rasterio.open(sys.argv[1]) as srast1:
    with rasterio.open(sys.argv[2]) as srast2:
        with rasterio.open(sys.argv[3]) as mapped:
            for y in range(image_height):
                srast1_scanline = srast1.read(1, window=Window(0, y, image_width, 1))
                srast2_scanline = srast2.read(1, window=Window(0, y, image_width, 1))
                mapped_scanline = mapped.read(1, window=Window(0, y, image_width, 1))

                for x in range(image_width):
                    val1 = int(srast1_scanline[0, x])
                    val2 = int(srast2_scanline[0, x])
                    map_val = int(mapped_scanline[0, x])
                    if (val1 < 0 or val2 < 0):
                        expected_map = -1
                    else:
                        expected_map = val1 * 49 + val2

                    if (expected_map != map_val):
                        raise ValueError("expected value: " + str(expected_map) + " real value: " + str(map_val) + " at index[" + str(x) + "][" + str(y) + "]")

                if (y % 100 == 0):
                    print(y)
                           
"""
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
"""
