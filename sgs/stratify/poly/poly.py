# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification using polygons
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

import tempfile

from sgs.utils import (
    SpatialRaster,
    SpatialVector,
)

from poly import poly_cpp

GIGABYTE = 1073741824
MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

def poly(
    raster: SpatialRaster,
    vector: SpatialVector,
    layer_name: str,
    attribute: str,
    features: list[str|list[str]],
    filename:str = '',
    driver_options: dict = None):
    """
    this function conducts stratification on a polygon by rasterizing a polygon
    layer, and using its values to determine stratifications.

    the layer_name parameter is the layer to be rasterized, and the attribute
    is the attribute within the layer to check. The features parameter specifies
    the feature values within the attribute, and which stratification they will
    be a part of.

    The features parameter is a list containing strings and lists of strings.
    The index within this list determines the stratification value. For example:
    
    features = ["low", "medium", "high"] 
        would result in 3 stratifications (0, 1, 2) where 'low' would correspond
        to stratification 0, medium to 1, and hight to 2

    features = ["low", ["medium", "high"]]
        would result in 2 stratifications (0, 1) where 'low' would correspond
        to stratification 0, and both medium and hight to 1

    Parameters
    --------------------
    rast : SpatialRaster
        raster data structure which will determine height, width, geotransform, and projection
    vector : SpatialVector
        the vector of polygons to stratify
    layer_name : str
        the layer in the vector to be stratified
    attribute : str
        the attribute in the layer to be stratified
    features : list[str|list[str]]
        the stratification values of each feature value, represented as the index in the list
    filename : str
        the output filename to write to, if desired

    Raises
    --------------------
    ValueError
        if the maximum strata value would result in an integer overflow error
    """

    cases = ""
    where_entries = []
    num_strata = len(features)

    if num_strata >= MAX_STRATA_VAL:
        raise ValueError("the number of features (and resulting max strata) will cause an overflow error because the max strata number is too large.")

    #generate query cases and where clause using features and attribute
    for i in range(len(features)):
        if type(features[i]) is not list:
            cases += "WHEN '{}' THEN {} ".format(str(features[i]), i)
            where_entries.append("{}='{}'".format(attribute, str(features[i])))
        else:
            for j in range(len(features[i])):
                cases += "WHEN '{}' THEN {} ".format(str(features[i][j]), i)
                where_entries.append("{}='{}'".format(attribute, str(features[i][j])))

    where_clause = " OR ".join(where_entries)

    #generate SQL query
    sql_query = f"""SELECT CASE {attribute} {cases}ELSE NULL END AS strata, {layer_name}.* FROM {layer_name} WHERE {where_clause}"""

    driver_options_str = {}
    if driver_options:
        for (key, val) in driver_options.items():
            if type(key) is not str:
                raise ValueError("the key for al key/value pairs in teh driver_options dict must be a string.")
            driver_options_str[key] = str(val)

    large_raster = raster.height * raster.width > GIGABYTE
    temp_dir = tempfile.mkdtemp()

    srast = SpatialRaster(poly_cpp(
        vector.cpp_vector,
        raster.cpp_raster,
        num_strata,
        layer_name,
        sql_query,
        filename,
        large_raster,
        temp_dir,
        driver_options_str
    ))

    #give srast ownership of it's own temp directory
    srast.have_temp_dir = True
    srast.temp_dir = temp_dir

    return srast

