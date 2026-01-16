# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification using polygons
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

##
# @defgroup user_poly poly
# @ingroup user_stratify

import os
import sys
import tempfile

from sgspy.utils import (
    SpatialRaster,
    SpatialVector,
)

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from _sgs import poly_cpp

GIGABYTE = 1073741824

##
# @ingroup user_poly
# This function conducts stratification on a vector dataset by rasterizing a polygon
# layer, and using its attribute values to determine stratifications.
# 
# the layer_name parameter is the layer to be rasterized, and the attribute
# is the attribute within the layer to check. The features parameter specifies
# the which feature value corresponds to which stratification.
# 
# The features parameter is a list containing strings and lists of strings.
# The index within this list determines the stratification value. For example:
# 
# features = ["low", "medium", "high"] @n
#     would result in 3 stratifications (0, 1, 2) where 'low' would correspond
#     to stratification 0, medium to 1, and hight to 2.
# 
# features = ["low", ["medium", "high"]] @n
#     would result in 2 stratifications (0, 1) where 'low' would correspond
#     to stratification 0, and both medium and high to stratification 1.
# 
# Examples
# --------------------
# rast = sgspy.SpatialRaster('rast.tif') @n
# vect = sgspy.SpatialVector('inventory_polygons.shp') @n
# srast = sgspy.stratify.poly(rast, vect, attribute='NUTRIENTS', layer_name='inventory_polygons', features=['poor', 'medium', 'rich'])
# 
# rast = sgspy.SpatialRaster('rast.tif') @n
# vect = sgspy.SpatialVector('inventory_polygons.shp') @n
# srast = sgspy.stratify.poly(rast, vect, attribute='NUTRIENTS', layer_name='inventory_polygons', 'features=['poor', ['medium', 'rich']], filename='nutrient_stratification.shp')
# 
# Parameters
# --------------------
# rast : SpatialRaster @n
#     raster data structure which will determine height, width, geotransform, and projection  @n @n
# vect : SpatialVector @n
#     the vector of polygons to stratify  @n @n
# layer_name : str @n
#     the layer in the vector to be stratified  @n @n
# attribute : str @n
#     the attribute in the layer to be stratified  @n @n
# features : list[str|list[str]] @n
#     the stratification values of each feature value, represented as the index in the list  @n @n
# filename : str @n
#     the output filename to write to, if desired  @n @n
# 
# Returns
# --------------------
# a SpatialRaster object containing the rasterized polygon.
def poly(
    rast: SpatialRaster,
    vect: SpatialVector,
    layer_name: str,
    attribute: str,
    features: list[str|list[str]],
    filename:str = '',
    driver_options: dict = None):

    MAX_STRATA_VAL = 2147483647 #maximum value stored within a 32-bit signed integer to ensure no overflow

    if type(rast) is not SpatialRaster:
        raise TypeError("'rast' parameter must be of type sgspy.SpatialRaster")

    if type(vect) is not SpatialVector:
        raise TypeError("'vect' parameter must be of type sgspy.SpatialVector")

    if type(layer_name) is not str:
        raise TypeError("'layer_name' parameter must be of type str.")

    if type(attribute) is not str:
        raise TypeError("'attribute' parameter must be of type str.")

    if type(features) is not list:
        raise TypeError("'features' parameter must be of type list.")

    if type(filename) is not str:
        raise TypeError("'filename' parameter must be of type str.")

    if driver_options is not None and type(driver_options) is not dict:
        raise TypeError("'driver_options' parameter, if givne, must be of type dict.")

    if rast.closed:
            raise RuntimeError("the C++ object which the rast object wraps has been cleaned up and closed.")

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

    large_raster = rast.height * rast.width > GIGABYTE
    
    #make temp directory which will be deleted if there is any problem when calling the cpp function
    temp_dir = tempfile.mkdtemp()
    rast.have_temp_dir = True
    rast.temp_dir = temp_dir

    srast = SpatialRaster(poly_cpp(
        vect.cpp_vector,
        rast.cpp_raster,
        num_strata,
        layer_name,
        sql_query,
        filename,
        large_raster,
        temp_dir,
        driver_options_str
    ))

    #now that it's created, give the cpp raster object ownership of the temporary directory
    rast.have_temp_dir = False
    srast.cpp_raster.set_temp_dir(temp_dir)
    srast.temp_dataset = filename == "" and large_raster
    srast.filename = filename

    return srast

