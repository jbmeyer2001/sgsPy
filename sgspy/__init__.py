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

def contains_all(path, files):
    missing = []
    for file in files:
        if not os.path.exists(os.path.join(path, file)):
            missing.append(file)

    return len(missing) == 0, missing

if platform.system() == 'Windows':
    #vendored gdal dll files and proj.db
    vendored_files = {"proj.db", "aec.dll", "charset-1.dll","freexl-1.dll","gdal.dll","geos.dll","geos_c.dll","geotiff.dll",
             "gif.dll","hdf5.dll","hdf5_cpp.dll","hdf5_hl.dll","hdf5_hl_cpp.dll","iconv-2.dll","jpeg62.dll","json-c.dll",
             "legacy.dll","Lerc.dll","libcrypto-3-x64.dll","libcurl.dll","libecpg.dll","libecpg_compat.dll","libexpat.dll",
             "liblzma.dll","libpgtypes.dll","libpng16.dll","libpq.dll","libsharpyuv.dll","libssl-3-x64.dll","libwebp.dll",
             "libwebpdecoder.dll","libwebpdemux.dll","libwebpmux.dll","libxml2.dll","lz4.dll","minizip.dll","netcdf.dll",
             "openjp2.dll","pcre2-16.dll","pcre2-32.dll","pcre2-8.dll","pcre2-posix.dll","proj_9.dll","qhull_r.dll",
             "spatialite.dll","sqlite3.dll","szip.dll","tiff.dll","turbojpeg.dll","uriparser.dll","zlib1.dll","zstd.dll"}
    
    root = os.path.dirname(os.path.realpath(__file__))

    """
    Check for the vendored binaries. These should be in the same folder as this file. If not, 
    it should be in site-packages (developer usage). If it is in neither, there is a problem.
    """
    found_all, missing = contains_all(root, vendored_files)
    if not found_all:
        root = os.path.join(list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0], "sgspy")
        found_all, _ = contains_all(root, vendored_files)

        #not in file path nor in expected environment location: throw error
        if not found_all:
            raise RuntimeError(f"{missing} not found. They should have been installed in the site-packages/sgspy directory of the current environment.")

        sys.path.append(root)
    os.environ["SGSPY_PROJDB_PATH"] = root

    #load all vendored dlls from correct place
    vendored_files.remove("proj.db")
    for file in vendored_files: ctypes.CDLL(os.path.join(root, file))

    """
    Check for external binaries. These are all a part of Python packages which have been installed,
    and so should be in a standard location for Python packages to put their binaries (depending on
    which environment manager is useb by the user).
    """
    external_dlls = ["onedal.3.dll", "tbb12.dll", "mkl_tbb_thread.2.dll", "mkl_core.2.dll", "onedal_core.3.dll", "onedal_thread.3.dll"]
    paths = [os.path.join(sys.prefix, "Library", "bin"),  
             os.path.join(site.USER_BASE, "Library", "bin"),
             os.path.join(sys.prefix, "DLLs")]

    for bin_path in paths:
        found_all, _ = contains_all(bin_path, external_dlls)
        if found_all:
            #load all external dlls from correct place
            for dll in external_dlls: ctypes.CDLL(os.path.join(bin_path, dll))
            break

    if not found_all:
        raise RuntimeError(f"""could not find a folder installed with all of the following dlls: {external_dlls}.
        checked the following paths: {paths}.
        First, ensure all of sgspy's Python dependencies are correctly installed.
        If they are, this is a bug and should be reported on https://github.com/jbmeyer2001/sgsPy/issues""")

else: #linux 
    root = os.path.dirname(os.path.realpath(__file__))

    """
    Check to ensure we can find proj.db
    """
    found_all, missing = contains_all(root, {"proj.db"})
    if not found_all:
        root = os.path.join(list(filter(lambda x : 'site-packages' in x, site.getsitepackages()))[0], "sgspy")
        found_all, _ = contains_all(root, {"proj.db"})

        #not in file path nor in expected environment location: throw error
        if not found_all:
            raise RuntimeError(f"{missing} not found. It should have been installed in the site-packages/sgspy directory of the current environment.")

        sys.path.append(root)
    os.environ["SGSPY_PROJDB_PATH"] = root
  
    #this library goes missing at runtime if we don't do this
    ctypes.CDLL(os.path.join(sys.prefix, 'lib', 'libtbb.so.12'), os.RTLD_GLOBAL | os.RTLD_NOW)

try:
    import _sgs
except ImportError as err:
    raise RuntimeError(f"""The following error has occured attempting to import _sgs: {err}.
    This is likely a bug, and should thus be reported on https://github.com/jbmeyer2001/sgsPy/issues""")

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
