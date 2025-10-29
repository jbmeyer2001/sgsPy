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
    This function conducts error-checking on user input parameters
    for systematic sampling, then calls a c++ function to conduct
    the sampling. The return values are plotted using matplotlib
    depending on the 'plot' parameter.

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
    plot : bool
        whether or not to plot the resulting samples
    filename : str
        the filename to write to or "" if not to write

    Raises
    --------------------
    ValueError
        if cellsize is less than or equal to 0
    ValueError
        if 'shape' parameter is not one of 'square', 'hexagon'
    ValueError
        if 'location' parameter is not on e of 'centers', 'corners', 'random'
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
