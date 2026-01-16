##
# @defgroup user User Documentation
# This is the documentation describing how to use the Python functions within the sgsPy
# package. For information on the underlying C++ implementations, see the developer
# docs.
#
# The first step in any processing using the sgsPy package will be to initialize in insance
# of either sgspy.SpatialRaster or sgspy.SpatialVector. These are the primary data inputs to
# all sgs functions, and information on their use can be found in the 'utils' section.
#
# The processing functions are split into three different categories: calculate, stratify,
# and sample. @n
# The calculate section contains various helpful functions to assist in sampling
# but are not necessarily a specific stratification or sampling function. Right now, 
# it only has 'pca' or principal component analysis. @n
# The stratify section has various stratification functions including stratification according 
# to user defined breaks 'breaks', stratification according to polygons 'poly', stratification
# along quantiles 'quantiles', and a method for mapping multiple existing stratificaiton outputs 'map'. @n
# The sample sections has various sampling functions including simple random sampling 'srs', stratified
# random sampling 'strat', systematic sampling 'systematic', and conditional latin hypercube sampling 'clhs'. @n

import os
import sys
import platform
import ctypes

if (platform.system() == 'Windows'):
    vendored_lib_path = os.path.join(sys.prefix, "sgspy")
    lib_path = os.path.join(sys.prefix, "Library", "bin")
    os.add_dll_directory(vendored_lib_path)
    os.add_dll_directory(lib_path)

    if vendored_lib_path not in os.environ['PATH']:
        os.environ['PATH'] = vendored_lib_path + os.pathsep + os.environ['PATH']

    if lib_path not in os.environ['PATH']:
        os.environ['PATH'] = lib_path + os.pathsep + os.environ['PATH']

else: #linux 
    #this library goes missing at runtime if we don't do this
    ctypes.CDLL(os.path.join(sys.prefix, 'lib', 'libtbb.so.12'), os.RTLD_GLOBAL | os.RTLD_NOW)
   
GIGABYTE = 1073741824

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
