# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification using polygons
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

from sgs.utils import (
    SpatialRaster,
    SpatialVector,
)

from poly import poly_cpp

def poly(
    vector: spatialVector,
    raster: spatialRaster,
    attribute: str,
    features: list[str|list[str]],
    filename:str = '',
    plot: bool = False):
    """

    """

    feature_lists = []

    for feature in features:
        if type(feature) is list:
            feature_lists.append(feature)
        else:
            feature_lists.append([feature])

    strat_rast = poly_cpp(
        vector.cpp_vector,
        raster.cpp_raster,
        attribute,
        feature_lists,
        filename
    )

    retval = SpatialRaster(strat_rast)

    if plot:
        retval.plot()

    return retval

