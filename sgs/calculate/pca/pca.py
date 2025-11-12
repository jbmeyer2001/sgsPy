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
    filename: str = '',
    display_info: bool = False,
    driver_options: dict = None
    ):
    """
    This functions conducts principal component analysis on the given
    raster.

    A number of output components must be provided as an integer. This integer
    must be less than or equal to the total number of bands in the input raster.
    This will be the number of bands in the output raster.
    A filename may be given to specify an output file location, otherwise
    a virtual file type will be used. The driver_options parameter is 
    used to specify creation options for a the output raster.
    See options for the Gtiff driver here: https://gdal.org/en/stable/drivers/raster/gtiff.html#creation-options

    Principal components are calculated across all raster bands, 
    along with mean and standard deviation of each raster band. The
    raster is both centered and scaled, then output values are calculated
    for each principal component.

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing input raster bands
    num_comp : int
        the number of components
    filename : str
        output filename or '' if there should not be an output file
    display_info : bool
        whether to display principal component eigenvalues/eigenvectors
    driver_options : dict
       the creation options as defined by GDAL which will be passed when creating output files

    Returns
    --------------------
    a SpatialRaster object containing principal component output bands.
    """
    breaks_dict = {}
    large_raster = False
    temp_folder = ""

    #ensure number of components is acceptabe
    if num_comp <= 0 or num_comp > len(rast.bands):
        msg = f"the number of components must be greater than zero and less than or equal to the total number of raster bands ({len(rast.bands)})."
        raise ValueError(msg)

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

    [pcomp, eigenvectors, eigenvalues] = pca_cpp(
        rast.cpp_raster,
        num_comp,
        large_raster,
        temp_dir,
        filename,
        driver_options_str
    )

    if display_info:
        print('eigenvectors:')
        print(eigenvectors)
        print()
        print('eigenvalues:')
        print(eigenvalues)
        print()

    pcomp_rast = SpatialRaster(pcomp)
    pcomp_rast.have_temp_dir = True
    pcomp_rast.temp_dir = temp_dir

    return pcomp_rast
