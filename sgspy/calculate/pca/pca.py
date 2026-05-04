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
import site
from sgspy.utils import SpatialRaster

from _sgs import pca_cpp

GIGABYTE = 1073741824

## 
# @ingroup user_pca
# This functions conducts principal component analysis on the given
# raster using the covariance matrix. Outputs are automatically
# scaled and shifted.
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
# pcomp, metadata = sgspy.calculate.pca(rast, 2, filename="pca.tif", return_metadata=True)
# eigenvectors, eigenvalues, means, stdevs = metadata
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
# return_metadata : bool @n
#     whether to return the eigenvectors, eigenvalues, means, and stdevs with the SpatialRaster @n @n
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
    return_metadata: bool = False,
    driver_options: dict = None
    ):
        
    if type(rast) is not SpatialRaster:
        print(type(rast))
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster.")

    if type(num_comp) is not int:
        raise TypeError("'num_comp' parameter must be of type int.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if type(return_metadata) is not bool:
        raise TypeError("'return_metadata' parameter must be of type bool.")

    if driver_options is not None and type(driver_options) is not dict: 
        raise TypeError("'driver_options' parameter, if given, must be of type dict.")

    if rast.closed:
        raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    breaks_dict = {}

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

    [pcomp, eigenvectors, eigenvalues, means, stdevs] = pca_cpp(
        rast.cpp_raster,
        num_comp,
        filename,
        driver_options_str
    )

    metadata = (eigenvectors, eigenvalues, means, stdevs)

    pcomp_rast = SpatialRaster(pcomp)

    if return_metadata:
        return pcomp_rast, metadata
    else:
        return pcomp_rast
