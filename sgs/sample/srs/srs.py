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
    SpatialRaster,
    SpatialVector,
    plot,
)

from srs import srs_cpp, srs_cpp_access

def srs(
    rast: SpatialRaster,
    num_samples: int,
    mindist: float = 0,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    plot: bool = False,
    filename: str = ''):
    """
    This function conducts simple random sampling on the raster given. 
    Sample points are randomly selected from data pixels (not be nodata).
    All sample points are at least mindist distance away from eachother.
    If unable to get the full number of sample points, a message is printed.

    An access vector of LineString or MultiLineString type can be provided.
    buff_outer specifies the buffer distance around the geometry which
    is allowed to be included in the sampling, buff_inner specifies the
    geometry which is not allowed to be included in the sampling. buff_outer
    must be larger than buff_inner. For a multi-layer vector, layer_name
    must be specified.

    Most of the calculation is done within the srs_cpp function which can
    be found in sgs/sample/srs/srs.cpp.

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to sample
    num_samples : int
        the target number of samples
    mindist : float
        the minimum distance each sample point must be from each other
    access (optional) : SpatialVector
        a vector specifying access network
    layer_name (optional) : str
        the layer within access that is to be used for sampling
    buff_inner (optional) : int | float
        buffer boundary specifying distance from access which CANNOT be sampled
    buff_outer (optional) : int | float
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
    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

    if (access):
        if layer_name is None:
            if len(access.layers) > 1:
                raise ValueError("if there are multiple layers in the access vector, layer_name must be defined.")
            layer_name = access.layers[0]

        if buff_inner is None:
            buff_inner = 0

        if buff_outer is None:
            raise ValueError("if an access vector is given, buff_outer must be defined.")

        if buff_inner >= buff_outer:
            raise ValueError("buff_outer must be greater than buff_inner")

    #call random sampling function
    if (access):
        [sample_coordinates, sample_points] = srs_cpp_access(
            rast.cpp_raster,
            num_samples,
            mindist,
            access.cpp_vector,
            layer_name,
            buff_inner,
            buff_outer,
            filename
        )
    else:
        [sample_coordinates, sample_points] = srs_cpp(rast.cpp_raster, num_samples, mindist, filename)
    
    if (len(sample_points)) < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, len(sample_points)))

    #plot new vector if requested
    #TODO do this in try/catch so that error won't cause
    #sampling to be thrown out
    if plot:
        fig, ax = plt.subplots()
        #TODO let user know which band is being printed
        rast.plot(ax, band=rast.bands[0])
        if access:
            access.plot('LineString', ax)
        ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
        plt.show()

    return sample_points
