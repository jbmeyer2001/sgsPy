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

from systematic import systematic_cpp

def systematic(
    rast: SpatialRaster,
    cellsize: float,
    shape: str = "square",
    location: str = "centers",
    existing: Optional[SpatialVector] = None,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    force: bool = False,
    plot: bool = False,
    filename: str = ""):
    """
    This function conducts systematic sampling within the extent of
    the raster given. The cellsize parameter specifies the grid size,
    shape specifies the grid type, and location specifies where in
    the grid a sample should fall into.

    shape can be one of 'square', and 'hexagon'.
    location can be one of 'corners', 'centers', 'random'.

    An access vector of LineString or MultiLineString type can be provided.
    buff_outer specifies the buffer distance around the geometry which is
    allowed to be included in the sampling, buff_inner specifies the geometry
    which is not allowed to be included in teh sampling. buff_outer must
    be larger than buff_inner. For a multi-layer vector, layer_name
    must be provided.

    A vector containing existing sample points can be provided. If this is
    the case then all of the points in the existing sample are automatically
    added and random samples are then chosen as required until num_samples 
    number of samples are chosen.

    If the force parameter is True, then the the samples are forced to 
    fall on an index which is NOT a no data value. This may result
    in some grids not being sampled.

    Parameters
    --------------------
    rast : SpatialRaster
        the raster to be sampled
    cellsize : float
        the size of the grid cells to be sampled
    shape : str
        the shape of the grid cells to be sampled
    location : str
        the location within the grid cell to be sampled
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
        whether or not to plot the resulting samples
    filename : str
        the filename to write to or "" if not to write
    """
    
    if cellsize <= 0:
        raise ValueError("cellsize must be greater than 0")

    if shape not in ["square", "hexagon"]:
        raise ValueError("shape parameter must be one of 'square', 'hexagon'")

    if location not in ["centers", "corners", "random"]:
        raise ValueError("location parameter must be one of 'centers', 'corners', 'random'")

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

    [samples, points, grid] = systematic_cpp(
        rast.cpp_raster,
        cellsize,
        shape,
        location,
        existing_vector,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        force,
        plot,
        filename
    )

    #plot new vector if requested
    if plot:
        fig, ax = plt.subplots()
        ax.set_xlim([rast.xmin, rast.xmax])
        ax.set_ylim([rast.ymin, rast.ymax])
        rast.plot(ax, band=rast.bands[0])
        title="samples on " + rast.bands[0]
        
        #plot grid
        for shape in grid:
            ax.plot(shape[0], shape[1], '-k')

        #plot sample points
        ax.plot(points[0], points[1], '.r')
        ax.set_title(label=title)
        plt.show()

    return SpatialVector(samples)
