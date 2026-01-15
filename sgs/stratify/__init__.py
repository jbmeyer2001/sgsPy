##
# @defgroup user_stratify stratify
# @ingroup user
#
# Documentation for the stratification functions. 

from . import (
    breaks,
    kmeans,
    poly,
    quantiles,
    map,
)

from .breaks import breaks
from .kmeans import kmeans
from .poly import poly
from .quantiles import quantiles
from .map import map

__all__ = [
    "breaks",
    "kmeans",
    "poly",
    "quantiles",
    "map",
]
