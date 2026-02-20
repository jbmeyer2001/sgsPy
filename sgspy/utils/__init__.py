##
# @defgroup user_utils utils
# @ingroup user
#
# Explanations of both the SpatialRaster and SpatialVector classes.

from . import (
    raster,
    vector,
)

from .raster import SpatialRaster
from .raster import StratRasterBandMetadata
from .vector import SpatialVector

__all__ = [
    "SpatialRaster",
    "StratRasterBandMetadata",
    "spatialVector",
]
