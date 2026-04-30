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
import site
import platform
import ctypes


potential_project_roots = [
    os.path.dirname(os.path.realpath(__file__)),
    os.path.join(list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0], "sgspy"),
]

for root in potential_project_roots:
    vendored_files_path = os.path.join(root, "vendored_files")

    if os.path.exists(vendored_files_path):
        if root not in os.environ["PATH"]: os.environ["PATH"] = root + os.pathsep + os.environ["PATH"]
        os.add_dll_directory(root)
        sys.path.append(root)
        
        if vendored_files_path not in os.environ["PATH"]: os.environ["PATH"] = vendored_files_path + os.pathsep + os.environ["PATH"]
        os.add_dll_directory(vendored_files_path)
        sys.path.append(vendored_files_path)
        os.environ["SGSPY_VENDORED_FILES_PATH"] = vendored_files_path
        
        break

if os.getenv("SGSPY_VENDORED_FILES_PATH") is None:
    raise RuntimeError("sgspy's vendored files path was unable to be found.")

if platform.system() == 'Windows':
    #ensure all dlls are able to be found
    for bin_path in [
        os.path.join(sys.prefix, "Library", "bin"), 
        os.path.join(sys.prefix, "DLLs"),
        os.path.join(site.USER_BASE, "Library", "bin")
    ]:    
        if os.path.exists(bin_path):
            os.add_dll_directory(bin_path)
            if bin_path not in os.environ["PATH"]: os.environ["PATH"] = bin_path + os.pathsep + os.environ["PATH"]

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
    StratRasterBandMetadata
)

from .calculate import (
    distribution,
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
