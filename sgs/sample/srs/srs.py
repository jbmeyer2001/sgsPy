# ******************************************************************************
#
#  Project: sgs
#  Purpose: simple random sampling (srs)
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import tempfile
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgs.utils import (
    SpatialRaster,
    SpatialVector,
    plot,
)

from srs import srs_cpp

def srs(
    rast: SpatialRaster,
    num_samples: int,
    mindist: float = 0,
    existing: Optional[SpatialVector] = None,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    plot: bool = False,
    filename: str = ''):
    """
    This function conducts simple random sampling on the raster given. 
    Sample points are randomly selected from data pixels (can't be nodata).
    All sample points are at least mindist distance away from eachother.
    If unable to get the full number of sample points, a message is printed.

    An access vector of LineString or MultiLineString type can be provided.
    buff_outer specifies the buffer distance around the geometry which
    is allowed to be included in the sampling, buff_inner specifies the
    geometry which is not allowed to be included in the sampling. buff_outer
    must be larger than buff_inner. For a multi-layer vector, layer_name
    must be specified.

    A vector containing existing sample points can be provided. If this is
    the case then all of the points in the existing sample are automatically
    added and random samples are chosen as required until num_samples number
    of samples are chosen.

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to sample
    num_samples : int
        the target number of samples
    mindist : float
        the minimum distance each sample point must be from each other
    existing (optional) : SpatialVector
        a vector specifying existing sample points
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
   """

    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

    if mindist is None:
        mindist = 0

    if mindist < 0:
        raise ValueError("mindist must be greater than or equal to 0")



    if (access):
        if layer_name is None:
            if len(access.layers) > 1:
                raise ValueError("if there are multiple layers in the access vector, layer_name parameter must be passed.")
            layer_name = access.layers[0]

        if layer_name not in access.layers:
            raise ValueError("layer specified by 'layer_name' does not exist in the access vector")

        if buff_inner is None or buff_inner < 0:
            buff_inner = 0

        if buff_outer is None or buff_outer < 0:
            raise ValueError("if an access vector is given, buff_outer must be a float greater than 0.")

        if buff_inner >= buff_outer:
            raise ValueError("buff_outer must be greater than buff_inner")

        access_vector = access.cpp_vector
    else:
        access_vector = None
        layer_name = ""
        buff_inner = -1
        buff_outer = -1

    if (existing):
        existing_vector = existing.cpp_vector
    else:
        existing_vector = None

    if not rast.have_temp_dir:
        rast.temp_dir = tempfile.mkdtemp()
        rast.have_temp_dir = True

    #call random sampling function
    [sample_coordinates, cpp_vector, num_points] = srs_cpp(
        rast.cpp_raster,
        num_samples,
        mindist,
        existing_vector,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        plot,
        rast.temp_dir,
        filename
    )
    
    if num_points < num_samples:
        print("unable to find the full {} samples within the given constraints. Sampled {} points.".format(num_samples, num_points))

    #plot new vector if requested
    if plot:
        try:
            fig, ax = plt.subplots()
            rast.plot(ax, band=rast.bands[0])
            title = "samples on " + rast.bands[0]
            
            if access:
                access.plot('LineString', ax)
                title += " with access"

            ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
            ax.set_title(label=title)
            plt.show()

        except Exception as e:
            print("unable to plot output: " + str(e))

    return SpatialVector(cpp_vector)
