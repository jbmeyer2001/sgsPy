# ******************************************************************************
#
#  Project: sgs
#  Purpose: balanced sampling
#  Author: Joseph Meyer
#  Date: July, 2025
#
# ******************************************************************************

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from sgs.utils import (
        SpatialRaster,
        SpatialVector,
        plot,
)

from balanced import (
        balanced_cpp, 
        balanced_strata_cpp,
        balanced_access_cpp,
        balanced_access_strata_cpp
)

def balanced(rast: SpatialRaster,
             num_samples: int,
             bands: Optional[list[str|int]] = None,
             algorithm: str = "lpm2_kdtree",
             srast: Optional[SpatialRaster] = None,
             srast_band: Optional[str|int] = None,
             prob: Optional[list[float]] = None,
             access: Optional[SpatialVector] = None,
             layer_name: Optional[str] = None,
             buff_inner: Optional[int | float] = None,
             buff_outer: Optional[int | float] = None,
             plot: bool = False,
             filename: str = ""):
    """
    This function conducts balanced sampling on the raster by calling a C++
    function, which makes use of the BalancedSampling R package C++ code.

    The lpm2_kdtree, lcube, and lcubestratified functionality are used, which 
    selects spatialy balanced samples. 

    For the lcubestratified functionality,
    an additional raster, srast, containing a stratification must be passed
    alongside srast_band, a string or integer indication of which band to
    use within the srast raster.

    the access parameter accepts a vector of Lines or Linestrings which
    represent road networks, rivers, etc. which impact the areas which
    are able to be sampled. The layer_name parameter represents the
    layer within the access vecter that should be used. The buff_inner
    parameter specifies the distance from the access netwrok which cannot be
    sampled within, the buff_outer paramter specifies the distance from
    the access network which must be sampled within.
    
    When the plot parameter is True, the resulting sample points will be
    plot on top of the raster and access vector if there is one.

    When the filename parameter is given, the output samples will be written
    to the given filename.

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure containing the raster to sample
    num_samples : int
        the target number of samples
    bands : Optional[list[str|int]]
        the bands in rast (either as a string or 0-indexed index value) to use for sampling
    algorithm : str
        indictation of balanced sampling method ('lmp2_kdtree', 'lcube', or 'lcubestratified')
    srast: Optional[SpatialRaster]
        must be given in algorithm is lcubestratified, the stratification raster
    srast_band : Optional[str|int]
        the band within srast which is the raster stratification
    prob : Optional[list[float]]
        the probabilities of each pixel for being in the final sample
    access : Optional[SpatialVector]
        vector of LineString or MultiLineString representing access network
    layer_name : Optional[str]
        the layer within the access vector to use
    buff_inner : Optional[int|float]
        buffer boundary specifying distance from access which CANNOT be sampled
    buff_outer : Optional[int|float]
        buffer boundary specifying distance from access which CAN be sampled
    plot : bool
        whether to plot the samples or not
    filename : str
        the filename to write to, or '' if file should not be written

    Raises
    --------------------
    ValueError
        if num_samples is less than 1
    ValueError
        if algorithm is not one of 'lpm2_kdtree', 'lcube', 'lcubestratified'
    ValueError
        if a band specified in teh bands list does not exist in the raster
    ValueError
        if lcubestratified is the algorithm but srast is not provided
    ValueError
        if lcubestratified is the algorithm and srast_band does not exist in srast
    ValueError
        if access vector is passed, and buff_outer is either not provided or less than or equal to 0
    ValueError
        if buff_inner is greater than buff_outer
    ValueError
        if access vector is passed, and layer_name is not in the access vector
    RuntimeError (C++)
        if allocating all the required memory for the BalancedSampling package is too large
    RuntimeError (C++)
        if the c++ code attempts to allocate data and is unable to
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
                    raise ValueError("band argument " + str(band) + "too large to be a zero-indexed band number")
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
    if prob is not None:
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
            [sample_coordinates, cpp_vector] = balanced_access_strata_cpp(
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
                plot,
                filename
            )
        else:
            [sample_coordinates, cpp_vector] = balanced_strata_cpp(
                rast.cpp_raster,
                num_samples,
                band_ints,
                srast.cpp_raster,
                srast_band,
                algorithm,
                prob,
                plot,
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
                plot,
                filename
            )
        else:
            [sample_coordinates, cpp_vector] = balanced_cpp(
                rast.cpp_raster,
                num_samples,
                band_ints,
                algorithm,
                prob,
                plot,
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
