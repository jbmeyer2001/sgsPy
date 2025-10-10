# ******************************************************************************
#
#  Project: sgs
#  Purpose: principal component analysis (pca)
#  Author: Joseph Meyer
#  Date: October, 2025
#
# ******************************************************************************

import tempfile
from sgs.utils import SpatialRaster
from pca import pca_cpp

GIGABYTE = 1073741824

def pca(
    rast: SpatialRaster,
    num_comp: int,
    plot: bool = False,
    filename: str = '',
    driver_options: dict = None
    ):
    """

    """
    breaks_dict = {}
    large_raster = False
    temp_folder = ""

    #ensure driver options keys are string, and convert driver options vals to string
    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise ValueError("the key for all key/value pairs in the driver_options dict must be a string.")
            driver_options_str[key] = str(val)

   #determine whether the raster should be categorized as 'large' and thus be processed in blocks
    raster_size_bytes = 0
    height = rast.height
    width = rast.width
    for i in range(len(rast.bands)):
        pixel_size = rast.cpp_raster.get_raster_band_type_size(i)
        band_size = height * width * pixel_size
        raster_size_bytes += band_size
        if band_size >= GIGABYTE:
            large_raster = True
            break

    large_raster = large_raster or (raster_size_bytes > GIGABYTE * 4)

    temp_dir = tempfile.mkdtemp()

    [pcomp, plot_results] = SpatialRaster(pca_cpp(
        rast.cpp_raster,
        num_comp,
        plot,
        large_raster,
        temp_dir,
        filename,
        driver_options_str
    ))
    
    pcomp_rast = SpatialRaster(pcomp)
    pcomp_rast.have_temp_dir = True
    pcomp_rast.temp_dir = temp_dir

    #TODO plot plot results???
    if plot:
        pass

    return pcomp_rast
    

