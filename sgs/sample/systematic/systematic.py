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
    
    if shape not in ["square", "hexagon", "triangle"]:
        raise RuntimeError("shape parameter must be one of 'square', 'hexagon', 'triangle'")

    if location not in ["centers", "corners", "random"]:
        raise RuntimeError("location parameter must be one of 'centers', 'corners', 'random'")

    [samples, points, grid] = systematic_cpp(
        rast.cpp_raster,
        cellsize,
        shape,
        location,
        plot,
        filename
    )

    #TODO delete this, it's for testing purposes
    fig, ax = plt.subplots()
    for shape in grid:
        ax.plot(shape[0], shape[1])
    plt.show()

    if plot:
        pass
        # plot points, grid on top of raster
        

