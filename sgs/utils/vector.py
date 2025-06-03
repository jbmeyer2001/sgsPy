import matplotlib.pyplot as plt

from vector import GDALVectorWrapper

class SpatialVector:
    """
    A wrapper class of a GDAL vector dataset.
    This class is primarily used under the hood, although it has a function
    for displaying vector information.

    Accessing vector info:
        
            vector metadata can be displayed using the info() function. All layers
            are displayed unless a specific layer is specified. The per-layer info
            includes: name, number of features, number of fields, geomtype, and bounds.

    Public Attributes:
    --------------------
    layer_names : list[str]
        a list of layer names
    
    Public Methods:
    --------------------
    info()
        takes an optional argument specify the band, and prints vector metadata to console
    """
    def __init__(self, image):
        """
        Constructing method for the SpatialVector class.

        Has one required parameter to specify a gdal dataset. The following
        attributes are populated:
        self.cpp_vector
        self.layer_names

        Parameters
        --------------------
        image: str
           specifies a path to a vector file

        Raises
        --------------------
        TypeError:
            if 'image' parameter is not of type str
        RuntimeError (from C++):
            if dataset is not initialized correctly 
        """
        if type(image) == str:
            self.cpp_vector = GDALVectorWrapper(image)
        else:
            raise TypeError(f"SpatialVector does not accept input of type {type(image)}")


        self.layers = self.cpp_vector.get_layer_names() 

    def print_info(self, name, layer_info):
        """
        prints layer information using the layer_info from self.cpp_vector.

        Parameters
        --------------------
        name : str
            str containing the layer name
        layer_info : dict
            dict containing 'feature_count', 'field_count', 'geometry_type', 'xmax', 'xmin', 'ymax', and 'ymin' items
        """
        print("{} layer info:".format(name))
        print("feature count: {}".format(layer_info['feature_count']))
        print("field count: {}".format(layer_info['field_count']))
        print("geometry type: {}".format(layer_info['geometry_type']))
        print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(
            layer_info['xmin'], 
            layer_info['xmax'], 
            layer_info['ymin'],
            layer_info['ymax']
        ))
        print()

    def info(self, layer=None):
        """
        calls self.print_info depending on layer parameter. If no layer is given,
        print all layers. A layer may be specified by either a str or an int.

        Parameters
        --------------------
        layer (optional) : str or int
            specifies the layer to print information on

        Raises
        --------------------
        TypeError:
            if the layer parameter is not of type None, str, or int
        """
        if layer is None:
            for layer in self.layers:
                self.print_info(layer, self.cpp_vector.get_layer_info(layer))
        elif type(layer) == str:
            self.print_info(layer, self.cpp_vector.get_layer_info(layer))
        elif type(layer) == int:
            self.print_info(layers[layer], self.cpp_vector.get_layer_info(layers[layer]))
        else:
              TypeError("layer parameter cannot be of type {}".format(type(layer)))
