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

from _sgs import clhs_cpp

def clhs(
    rast: SpatialRaster,
    num_samples: int,
    iterations: int = 10000,
    access: Optional[SpatialVector] = None,
    layer_name: Optional[str] = None,
    buff_inner: Optional[int | float] = None,
    buff_outer: Optional[int | float] = None,
    plot: bool = False,
    filename: str = ''):
    """

    """
    if rast.closed:
            raise RuntimeError("the C++ object which the raster object wraps has been cleaned up and closed.")

    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

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

    temp_dir = rast.cpp_raster.get_temp_dir()
    if temp_dir == "":
        temp_dir = tempfile.mkdtemp()
        rast.cpp_raster.set_temp_dir(temp_dir)

    [sample_coordinates, cpp_vector] = clhs_cpp(
        rast.cpp_raster,
        num_samples,
        iterations,
        access_vector,
        layer_name,
        buff_inner,
        buff_outer,
        plot,
        temp_dir,
        filename
    )

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
