# ******************************************************************************
#
#  Project: sgs
#  Purpose: principal component analysis (pca)
#  Author: Joseph Meyer
#  Date: October, 2025
#
# ******************************************************************************

##
# @defgroup user_pca pca
# @ingroup user_calculate

import os
import sys
import tempfile
from sgspy.utils import SpatialRaster

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from _sgs import pca_cpp

GIGABYTE = 1073741824

## 
# @ingroup user_pca
# This functions conducts principal component analysis on the given
# raster.
# 
# A number of output components must be provided as an integer. This integer
# must be less than or equal to the total number of bands in the input raster,
# and will be the number of bands in the output raster.
# A filename may be given to specify an output file location, otherwise
# a virtual file type will be used. The driver_options parameter is 
# used to specify creation options for a the output raster.
# See options for the Gtiff driver here: https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options
# 
# Principal components are calculated across all raster bands, 
# along with mean and standard deviation of each raster band. The
# raster is both centered and scaled, then output values are calculated
# for each principal component.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster("raster.tif") @n
# pcomp = sgspy.calculate.pca(rast, 3) 
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# pcomp = sgspy.calculate.pca(rast, 2, filename="pca.tif", display_info=True)
# 
# rast = sgspy.SpatialRaster("raster.tif") @n
# pcomp = sgspy.calculate.pca(rast, 1, filename="pca.tif", driver_options={"COMPRESS": "LZW"}) 
# 
# Parameters
# --------------------
# rast : SpatialRaster @n
#     raster data structure containing input raster bands @n @n
# num_comp : int @n
#     the number of components @n @n
# filename : str @n
#     output filename or '' if there should not be an output file @n @n
# display_info : bool @n
#     whether to display principal component eigenvalues/eigenvectors @n @n
# driver_options : dict @n
#    the creation options as defined by GDAL which will be passed when creating output files @n @n
# 
# Returns
# --------------------
# a SpatialRaster object containing principal component bands
def pca(
    rast: SpatialRaster,
    num_comp: int,
    filename: str = '',
    display_info: bool = False,
    driver_options: dict = None
    ):
        
    if type(rast) is not SpatialRaster:
        print(type(rast))
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster.")

    if type(num_comp) is not int:
        raise TypeError("'num_comp' parameter must be of type int.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if type(display_info) is not bool:
        raise TypeError("'display_info' parameter must be of type bool.")

    if driver_options is not None and type(driver_options) is not dict: 
        raise TypeError("'driver_options' parameter, if given, must be of type dict.")

    if rast.closed:
        raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    breaks_dict = {}
    large_raster = False
    temp_folder = ""

    #ensure number of components is acceptabe
    if num_comp <= 0 or num_comp > len(rast.bands):
        msg = f"the number of components must be greater than zero and less than or equal to the total number of raster bands ({len(rast.bands)})."
        raise ValueError(msg)

    #ensure driver options keys are string, and convert driver options vals to string
    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise TypeError("the key for all key/value pairs in the driver_options dict must be a string.")
            driver_options_str[key] = str(val)

   #determine whether the raster should be categorized as 'large' and thus be processed in blocks
    raster_size_bytes = 0
    height = rast.height
    width = rast.width
    for i in range(len(rast.bands)):
        pixel_size = rast.cpp_raster.get_raster_band_type_size(i)
        band_size = height * width * pixel_size
        raster_size_bytes += band_size
        if band_size >= GIGABYTE:
            large_raster = True
            break

    large_raster = large_raster or (raster_size_bytes > GIGABYTE * 4)

    temp_dir = tempfile.mkdtemp()

    [pcomp, eigenvectors, eigenvalues] = pca_cpp(
        rast.cpp_raster,
        num_comp,
        large_raster,
        temp_dir,
        filename,
        driver_options_str
    )

    if display_info:
        print('eigenvectors:')
        print(eigenvectors)
        print()
        print('eigenvalues:')
        print(eigenvalues)
        print()

    pcomp_rast = SpatialRaster(pcomp)
    pcomp_rast.have_temp_dir = True
    pcomp_rast.temp_dir = temp_dir
    pcomp_rast.temp_dataset = filename == "" and large_raster
    pcomp_rast.filename = filename

    return pcomp_rast
