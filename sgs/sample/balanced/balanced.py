from typing import Optional

from sgs.utils import (
        access,
        SpatialRaster,
        SpatialVector,
        plot,
        write,
)

from balanced import test_func

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

    print("calling test_func on raster")
    test_func(raster)
    print("called it")
    return

    if algorithm not in ["lpm2_kdtree", "lcube", "lcubestratified"]:
        raise ValueError("algorithm parameter must specify one of: 'lpm2_kdtree', 'lcube', 'lcubestratified'.")

    if access:
        #TODO add when access has been implemented
        sgs.utils.access()

    num_pixels = raster.width * raster.height
    if not prob:
        prob = [num_samples / num_pixels] * num_pixels

    if algorithm == "lpm2_kdtree":
        samples = None
        #call C++ bound function

    if algorithm == "lcube":
        samples = None
        #call C++ bound function

    if algorithm == "lcubestratified":
        samples = None
        #call C++ bound function
        #this one will likely require some checking to ensure there is a 'strata' layer in the raster
    
    #TODO: convert coordinates to spatial points

    if plot:
        #TODO add when plot has been implemented
        sgs.utils.plot()

    if filename:
        #TODO add when write has been implemented
        sgs.utils.write(filename, overwrite)

    print(__file__)
    raise NotImplementedError
