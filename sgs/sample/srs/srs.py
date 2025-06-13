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
import matplotlib.pyplot as plt

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
    [sample_coordinates, sample_points] = srs_cpp(rast.cpp_raster, mindist, num_samples, filename)
    
    if (len(sample_points)) < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, len(sample_points)))

    if plot:
        fig, ax = plt.subplots()
        rast.plot(ax)
        if access is not None:
            access.plot('LineString', ax)
        ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
        plt.show()

    if filename:
        #TODO add when write has been implemented
        write(filename, overwrite)

    return sample_points
