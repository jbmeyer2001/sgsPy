# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
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
    This function conducts simple random sampling on the raster given. 
    Sample points are randomly selected from data pixels (can't be nodata).
    All sample points are at least mindist distance away from eachother.
    If unable to get the full number of sample points, a message is printed.

    Most of the calculation is done within the srs_cpp function which can
    be found in sgs/sample/srs/srs.cpp.

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to sample
    num_samples : int
        the target number of samples
    access (optional) : SpatialVector
        a vector specifying access information
    mindist : float
        the minimum distance each sample point must be from each other
    buf_inner (optional) : int | float
        buffer boundary specifying distance from access which CANNOT be sampled
    buf_outer (optional) : int | float
        buffer boundary specifying distance from access which CAN be sampled
    plot : bool
        whether to plot the samples or not
    filename : str
        the filename to write to, or '' if file should not be written

    Raises
    --------------------
        ValueError
            if num_samples is less than 1
        RuntimeError (from C++)
            if the number of samples is greater than the number of data pixels in the image
        RuntimeError (from C++)
            if there is an issue reading a band from the raster
        RuntimeError (from C++)
            type errors for not valid/accepted types -- if this error shows up from this function, it means there's a bug
    """
    if access is not None:
        pass
        #call access function TODO

    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

    #call random sampling function
    [sample_coordinates, sample_points] = srs_cpp(rast.cpp_raster, mindist, num_samples, filename)
    
    if (len(sample_points)) < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, len(sample_points)))

    #plot new vector if requested
    if plot:
        fig, ax = plt.subplots()
        rast.plot(ax)
        if access is not None:
            access.plot('LineString', ax)
        ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
        plt.show()

    return sample_points
