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

    """
    
    if cellsize <= 0:
        raise ValueError("cellsize must be greater than 0")

    if shape not in ["square", "hexagon", "triangle"]:
        raise ValueError("shape parameter must be one of 'square', 'hexagon', 'triangle'")

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
