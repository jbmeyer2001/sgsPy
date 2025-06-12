from typing import Optional

import numpy as np

from sgs.utils import (
    access,
    SpatialRaster,
    SpatialVector,
    plot,
    write,
)


def srs(
    rast: SpatialRaster,
    num_samples: int,
    access: Optional[SpatialVector] = None,
    mindist: Optional[float] = None,
    buf_inner: Optional[int | float] = None,
    buf_outer: Optional[int | float] = None,
    plot: bool = false,
    filename: Optional[str] = None,
    overwrite: bool = False):
    """

    """

    if access is not None:
        pass
        #call access function TODO

    #call cpp srs function
    samples = srs_cpp()

    if plot:
        pass
        #call plot function

    if filename:
        #TODO add when write has been implemented
        write(filename, overwrite)

    return samples
