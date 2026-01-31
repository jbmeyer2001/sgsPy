##
# @defgroup user_calculate calculate
# @ingroup user
#
# documentation of additional calculation functions for sgsPy. At the moment just principal component analysis.

from . import (
    distribution,
    pca,
    representation,
)

from .distribution import distribution
from .pca import pca
from .representation import representation

__all__ = [
    "distribution",
    "pca",
    "representation",
]
