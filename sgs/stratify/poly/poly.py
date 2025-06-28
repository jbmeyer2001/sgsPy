# ******************************************************************************
#
#  Project: sgs
#  Purpose: stratification using polygons
#  Author: Joseph Meyer
#  Date: June, 2025
#
# ******************************************************************************

from sgs.utils import (
    SpatialRaster,
    SpatialVector,
)

from poly import poly_cpp

def poly(
    raster: SpatialRaster,
    vector: SpatialVector,
    attribute: str,
    layer_name: str,
    features: list[str|list[str]],
    filename:str = '',
    plot: bool = False):
    """

    note that I'm using the GDAL rasterize sql functionality because it
    can all be condensed into one, somewhat concise query. No columns
    will need to be added to the dataset, and no dataset parameters
    will have to be modified

    """

    cases = ""
    where_entries = []
    #generate query cases and where clause using features and attribute
    for i in range(len(features)):
        if type(features[i]) is str:
            cases += "WHEN '{}' THEN {} ".format(features[i], i)
            where_entries.append("{}='{}'".format(attribute, features[i]))
        else:
            for j in range(len(features[i])):
                cases += "WHEN '{}' THEN {} ".format(features[i][j], i)
                where_entries += "{}='{}'".format(attribute, features[i][j])

    where_clause = " OR ".join(where_entries)

    #generate SQL query
    sql_query = f"""SELECT CASE {attribute} {cases}ELSE NULL END AS strata, {layer_name}.* FROM {layer_name} WHERE {where_clause}"""

    strat_rast = poly_cpp(
        vector.cpp_vector,
        raster.cpp_raster,
        sql_query,
        filename
    )

    retval = SpatialRaster(strat_rast)

    if plot:
        retval.plot()

    return retval

