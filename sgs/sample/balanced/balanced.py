from typing import Optional

import numpy as np

from sgs.utils import (
        SpatialRaster,
        SpatialVector,
        plot,
)
from balanced import (
        lcube_cpp, 
        lcube_stratified_cpp, 
        hlpm2_cpp
)
from balanced import (
        lcube_cpp, 
        lcube_stratified_cpp, 
        hlpm2_cpp
)

#optional srast for lcube_stratified ?????
def balanced(rast: SpatialRaster,
             num_samples: int,
             bands = Optional[list[str]] = None;
             algorithm: str = "lpm2_kdtree",
             prob: Optional[list[float]] = None,
             access: Optional[SpatialVector] = None,
             layer_name: Optional[str] = None
             buf_inner: Optional[int | float] = None,
             buf_outer: Optional[int | float] = None,
             plot: bool = False,
             filename: Optional[str] = None,
             overwrite: bool = False):
    """
    Balanced sampling using #### package ...
    (add documentation)

    something about algorithm constraints

    something about access
    """

    if algorithm not in ["lpm2_kdtree", "lcube", "lcubestratified"]:
        raise ValueError("algorithm parameter must specify one of: 'lpm2_kdtree', 'lcube', 'lcubestratified'.")

    if prob:
        prob = np.ascontiguousarray(
            prob, 
            dtype=np.float64
        )
    else:
        prob = np.ascontiguousarray([], dtype=np.float64)

    

    if algorithm == "lpm2_kdtree":
        samples = hlpm2_cpp(rast.cpp_raster, prob.data)
    elif algorithm == "lcube":
        samples = lcube_cpp(rast.cpp_raster, prob.data)
    else:
        #change this as the band may not be named 'strata' exactly
        if 'strata' not in raster.bands:
            raise ValueError("raster must have a band 'strata' if using lcubestratified method.")
        samples = lcube_stratified_cpp(rast.cpp_raster, prob.data)

    if plot:
        #TODO add when plot has been implemented
        sgs.utils.plot()

    return samples
