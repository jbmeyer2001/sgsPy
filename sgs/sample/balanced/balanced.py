from typing import Optional
from osgeo import gdal

from sgs.utils import (
        access,
        SpatialRaster,
        SpatialVector,
        plot,
        write,
)

def balanced(raster: SpatialRaster,
             num_samples: int,
             algorithm: str = "lpm2_kdtree",
             prob: Optional[list[float]] = None,
             access: Optional[SpatialVector] = None,
             buf_inner: Optional[int | float] = None,
             buf_outer: Optional[int | float] = None,
             plot: bool = False,
             filename: Optional[str] = None,
             overwrite: bool = False):
    """
    Balanced sampling using #### package ...
    (add documentation)

    something about algorithm constraints

    something about access
    """
    
    vmemaddr = raster.dataset.GetVirtualMem(
        eRWFlag=gdal.GA_Update,
        nXOff=0,
        nYOff=0,
        nXSize=raster.width,
        nYSize=raster.height,
        nBufXSize=raster.width,
        nBufYSize=raster.height,
        eBufType=raster.dataset.GetRasterBand(1).DataType,
        band_list=list(range(1, raster.layers + 1)),
        bIsBandSequential=True,
        nCacheSize=10 * 1024 * 1024, #TODO this is from documentation, likely will be changed
        nPageSizeHint=0,
        options={}
    ).GetAddr()
    
    match raster.data_type:
        case gdal.GDT_Int8:
            cpp_type = 'int8_t'
        case gdal.GDT_UInt16:
            cpp_type = 'uint16_t'
        case gdal.GDT_Int16:
            cpp_type = 'int16_t'
        case gdal.GDT_UInt32:
            cpp_type = 'uint32_t'
        case gdal.GDT_Int32:
            cpp_type = 'int32_t'
        case gdal.GDT_Float32:
            cpp_type = 'float32_t'
        case gdal.GDT_Float64:
            cpp_type = 'float64_t'
        case gdal.GDT_UInt64:
            cpp_type = 'uint64_t'
        case gdal.GDT_Int64:
            cpp_type = 'int64_t'
        case _:
            raise TypeError("""
            Raster pixel type is not one of the acceptable types. 
            The acceptable types include: 
                GDT_Int8
                GDT_UInt16
                GDT_Int16
                GDT_UInt32
                GDT_Int32
                GDT_Float32
                GDT_Float64
                GDT_UInt64
                GDT_Int64
            """)

    print('memory address: ' + str(vmemaddr))
    print('cpp type: ' + str(cpp_type))
    return

    if algorithm not in ["lpm2_kdtree", "lcube", "lcubestratified"]:
        raise ValueError("algorithm parameter must specify one of: 'lpm2_kdtree', 'lcube', 'lcubestratified'.")

    if access:
        #TODO add when access has been implemented
        sgs.utils.access()

    num_pixels = raster.width * raster.height
    if not prob:
        prob = [num_samples / num_pixels] * num_pixels

    if algorithm == "lpm2_kdtree":
        samples = None
        #call C++ bound function

    if algorithm == "lcube":
        samples = None
        #call C++ bound function

    if algorithm == "lcubestratified":
        samples = None
        #call C++ bound function
        #this one will likely require some checking to ensure there is a 'strata' layer in the raster
    
    #TODO: convert coordinates to spatial points

    if plot:
        #TODO add when plot has been implemented
        sgs.utils.plot()

    if filename:
        #TODO add when write has been implemented
        sgs.utils.write(filename, overwrite)

    print(__file__)
    raise NotImplementedError
