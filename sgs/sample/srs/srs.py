# ******************************************************************************
#
#  Project: sgs
#  Purpose: GDALDataset wrapper for raster operations
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

from typing import Optional

import numpy as np

from sgs.utils import (
    access,
    SpatialRaster,
    SpatialVector,
    plot,
    write,
)

from srs import srs_cpp

def srs(
    rast: SpatialRaster,
    num_samples: int,
    access: Optional[SpatialVector] = None,
    mindist: float = 0,
    buf_inner: Optional[int | float] = None,
    buf_outer: Optional[int | float] = None,
    plot: bool = False,
    filename: str = ''):
    """

    """
    if access is not None:
        pass
        #call access function TODO

    #call random sampling function
    [sample_coordinates, sample_wkt] = srs_cpp(rast.cpp_raster, mindist, num_samples, filename)
    print(sample_coordinates)
    print(sample_wkt)
    if plot:
        pass
        #call plot function

    if filename:
        #TODO add when write has been implemented
        write(filename, overwrite)

    return samples
