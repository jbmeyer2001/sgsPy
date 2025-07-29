from typing import Optional

import numpy as np

from sgs.utils import (
        SpatialRaster,
        SpatialVector,
        plot,
)
from balanced import (
        lcube_cpp, 
        lcube_stratified_cpp, 
        hlpm2_cpp
)
from balanced import (
        balanced_cpp, 
        balanced_strat_cpp,
        balanced_access_cpp,
        balanced_access_strat_cpp
)

#optional srast for lcube_stratified ?????
def balanced(rast: SpatialRaster,
             num_samples: int,
             bands = Optional[list[str|int]] = None;
             algorithm: str = "lpm2_kdtree",
             srast: Optional[SpatialRaster] = None,
             srast_band: Optional[str|int] = None,
             prob: Optional[list[float]] = None,
             access: Optional[SpatialVector] = None,
             layer_name: Optional[str] = None
             buf_inner: Optional[int | float] = None,
             buf_outer: Optional[int | float] = None,
             plot: bool = False,
             filename: Optional[str] = None,
             overwrite: bool = False):
    """
    Balanced sampling using #### package ...
    (add documentation)

    something about zero indexed bands

    something about algorithm constraints

    something about access
    """
    if num_samples < 1:
        raise ValueError("num_samples must be greater than 0")

    #algorithm argument error checking
    if algorithm not in ["lpm2_kdtree", "lcube", "lcubestratified"]:
        raise ValueError("algorithm parameter must specify one of: 'lpm2_kdtree', 'lcube', 'lcubestratified'")

    #band argument error checking
    band_ints = []
    if bands is not None:
        for band in bands:
            if type(band) == int:
                if band > len(rast.bands) - 1:
                    raise ValueError("band argument " + str(band) + "too large to be a zero-indexed band number"
                else:
                    band_ints.append(band)
            else: #type(bands) == str
                if band not in rast.bands:
                    raise ValueError("band argument " + band + " is not one of the raster bands")
                else:
                    band_ints.append(rast.band_name_dict[band])
    else:
        band_ints = [i for i in range(len(rast.bands))]

    #prob argument modification
    if prob:
        prob = np.ascontiguousarray(prob, dtype=np.float64)
    else:
        prob = np.ascontiguousarray([], dtype=np.float64)

    #srast and srast_band argument checking
    if algorithm == "lcubestratified":
        if srast is None:
            raise ValueError("srast argument must be provided if the algorithm is lcubestratified")
        elif srast_band is None:
            raise ValueError("srast_band argument must be provided if the algorithm is lcubestratified")
        elif type(srast_band) == int and srast_band > len(srast.bands) - 1:
            raise ValueError("srast_band " + str(srast_band) + " is not a valid band index in srast")
        elif type(srast_band) == str and srast_band not in srast.bands:
            raise ValueError("srast_band " + srast_band + " is not one of the band names in srast")

        #convert to int if string
        if type(srast_band) == str:
            srast_band = srast.band_name_dict[srast_band]

    #access argument checking
    if (access):
        if layer_name is None:
            if len(access.layers) > 1:
                raise ValueError("if there are multiple layers in the access vector, layer_name parameter must be passed.")
            layer_name = access.layers[0]

        if layer_name not in access.layers:
            raise ValueError("layer specified by 'layer_name' does not exist in the access vector")

        if buff_inner is None or buff_inner < 0:
            buff_inner = 0

        if buff_outer is None or buff_outer < 0:
            raise ValueError("if an access vector is given, buff_outer must be a float greater than 0.")

        if buff_inner >= buff_outer:
            raise ValueError("buff_outer must be greater than buff_inner")

    if algorithm == "lcubestratified":
        if access:
            [sample_coordinates, cpp_vector] = balanced_access_strat_cpp(
                rast.cpp_raster,
                num_samples,
                band_ints,
                srast.cpp_raster,
                srast_band,
                access.cpp_vector,
                layer_name,
                buff_inner,
                buff_outer,
                algorithm,
                prob,
                filename
            )
        else:
            [sample_coordinates, cpp_vector] = balanced_strat_cpp(
                rast.cpp_raster,
                num_samples,
                band_ints,
                srast.cpp_raster,
                srast_band,
                algorithm,
                prob,
                filename
            )
    else:
        if access:
            [sample_coordinates, cpp_vector] = balanced_access_cpp(
                rast.cpp_raster,
                num_samples,
                band_ints,
                access.cpp_vector,
                layer_name,
                buff_inner,
                buff_outer,
                algorithm,
                prob,
                filename
            )
        else:
            [sample_coordinates, cpp_vector] = balanced_cpp(
                rast.cpp_raster,
                num_samples,
                band_ints,
                algorithm,
                prob,
                filename,
            )

    #plot new vector if requested
    if plot:
        try:
            fig, ax = plt.subplots()
            rast.plot(ax, band=rast.bands[band_ints[0]])
            title = "samples in " + rast.bands[band_ints[0]]

            if access:
                access.plot('LineString', ax)
                title += " with access"

            ax.plot(sample_coordinates[0], sample_coordinates[1], '.r')
            ax.set_title(label=title)
            plt.show()

        except Exception as e:
            print("unable to plot output: " + str(e))

    return SpatialVector(cpp_vector)
