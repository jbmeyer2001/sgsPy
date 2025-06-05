from typing import Optional

from sgs.utils import (
        access,
        SpatialRaster,
        SpatialVector,
        plot,
        write,
)

from balanced import lcube_cpp, lcube_stratified_cpp, hlpm2_cpp

def balanced(raster: SpatialRaster,
             num_samples: int,
             algorithm: str = "lpm2_kdtree",
             prob: Optional[list[float]] = None,
             access: Optional[SpatialVector] = None,
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

    if access:
        access_vector = access.cpp_vector
    else:
        access_vector = None

    if algorithm == "lpm2_kdtree":
        samples = hlpm2_cpp(raster.cpp_raster, access_vector, prob)

    if algorithm == "lcube":
        samples = lcube_cpp(raster.cpp_raster, access_vector, prob)

    if algorithm == "lcubestratified":
        if 'strata' not in raster.bands:
            raise ValueError("raster must have a band 'strata'")
        samples = lcube_stratified_cpp(raster.cpp_raster, access_vector, prob)
    
    #TODO: convert coordinates to spatial points

    if plot:
        #TODO add when plot has been implemented
        sgs.utils.plot()

    if filename:
        #TODO add when write has been implemented
        sgs.utils.write(filename, overwrite)

    print(__file__)
    raise NotImplementedError
