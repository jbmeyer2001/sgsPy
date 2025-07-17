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
    plot,
)

from systematic import systematic_cpp

def systematic(
    rast: SpatialRaster,
    cellsize: float,
    shape: str = "square",
    location: str = "centers",
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

    [samples, points, grid] = systematic_cpp(
        rast.cpp_raster,
        cellsize,
        shape,
        location,
        plot,
        filename
    )

    if plot:
        fig, ax = plt.subplots()
        ax.set_xlim([rast.xmin, rast.xmax])
        ax.set_ylim([rast.ymin, rast.ymax])

        #TODO let user know which band is being printed
        #plot raster
        rast.plot(ax, band=rast.bands[0])
        
        #plot grid
        for shape in grid:
            ax.plot(shape[0], shape[1], '-k')

        #plot sample points
        ax.plot(points[0], points[1], '.r')
        plt.show()

    return samples
