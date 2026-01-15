##
# @defgroup user_calculate calculate
# @ingroup user
#
# documentation of additional calculation functions for sgsPy. At the moment just principal component analysis.

from . import (
    pca,
    representation,
)

from .pca import pca
from .representation import representation

__all__ = [
    "pca",
    "representation",
]
