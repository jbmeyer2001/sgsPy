import os
import sys
import platform
import ctypes

if (platform.system() == 'Windows'):
    bin_path = os.path.join(sys.prefix, "sgs")
    os.add_dll_directory(bin_path)

    if bin_path not in os.environ.get('PROJ_LIB', ''):
        os.environ['PROJ_LIB'] = bin_path

    if bin_path not in os.environ['PATH']:
        os.environ['PATH'] = bin_path + os.pathsep + os.environ['PATH']

else: #linux 
    #this library goes missing at runtime if we don't do this
    ctypes.CDLL(os.path.join(sys.prefix, 'lib', 'libtbb.so.12'), os.RTLD_GLOBAL | os.RTLD_NOW)
    
from . import utils
from . import calculate
from . import sample
from . import stratify

from .utils import (
    SpatialRaster,
    SpatialVector,
)

from .calculate import (
    pca,
    representation,
)

from .sample import (
    ahels,
    clhs,
    nc,
    srs,
    strat,
    systematic,
)

from .stratify import (
    breaks,
    kmeans,
    poly,
    quantiles,
    map,
)

__all__ = list(
    set(utils.__all__) |
    set(calculate.__all__) |
    set(sample.__all__) |
    set(stratify.__all__)
)
