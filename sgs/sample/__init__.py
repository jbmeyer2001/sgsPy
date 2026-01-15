## 
# @defgroup user_sample sample
# @ingroup user
#
# Documentation for the sampling functions.

from . import (
    ahels,
    clhs,
    nc,
    srs,
    strat,
    systematic,
)

from .ahels import ahels
from .clhs import clhs
from .nc import nc
from .srs import srs
from .strat import strat
from .systematic import systematic

__all__ = [
    "ahels",
    "clhs",
    "nc",
    "srs",
    "strat",
    "systematic",
]
