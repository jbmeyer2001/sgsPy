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

from systematic import systematic_cpp, systematic_cpp_access

def systematic(
    rast: SpatialRaster,
    cellsize: float,
    shape: str = "square",
    location: str = "centers",
    plot: bool = false,
    filename: str = ""):
    """

    """
    
    if shape not in ["square", "hexagon", "triangle"]:
        raise RuntimeError("shape parameter must be one of 'square', 'hexagon', 'triangle'")

    if location not in ["centers", "corners", "random"]:
        raise RuntimeError("location parameter must be one of 'centers', 'corners', 'random'")

    [samples, points, grid] = systematic_cpp(
        rast.cpp_rast,
        cellsize,
        shape,
        location,
        plot,
        filename
    )

    if plot:
        # plot points, grid on top of raster
        

