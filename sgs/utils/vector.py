from osgeo import gdal, ogr
import matplotlib.pyplot as plt

class SpatialVector:
    """
    A wrapper class of a GDAL vector dataset.
    This class is primarily used under the hood, althoufh it has a function
    for displaying vector information.

    Accessing vector info:
        
            vector metadata can be displayed using the info() function. All layers
            are displayed unless a specific layer is specified. The per-layer info
            includes: name, number of features, number of fields, geomtypes, and bounds.

    Public Attributes:
    --------------------
    layer_count : int
        the number of layers in the vector image
    layers_info : list[dict]
        a list of layers with dicts containing name, features, fields, geomtypes, and extent info
    
    Public Methods:
    --------------------
    info()
        takes no arguments, prints vector metadata to console
    """
    def __init__(self, image):
        """
        Constructing method for the SpatialVector class.

        Has one required parameter to specify a gdal dataset. The following
        attributes are populated using the given dataset:
        self.dataset
        self.layer_count
        self.layers
        self.layers_info
        self.layer_name_to_index

        Parameters
        --------------------
        image: str OR gdal.Dataset
           specifies either a gdal dataset, or a path to a gdal dataset 

        Raises
        --------------------
        TypeError:
            if 'image' parameter is not of type str or gdal.Dataset
        ValueError:
            if dataset is not loaded
        """
        if type(image) == str:
            self.dataset = gdal.OpenEx(image, gdal.OF_VECTOR)
        elif type(image) == gdal.Dataset:
            self.dataset = image
        else:
            raise TypeError(f"SpatialVector does not accept input of type {type(image)}")

        if not self.dataset:
            raise ValueError("dataset must exist")

        self.layer_count = self.dataset.GetLayerCount()
        self.layers = []
        self.layers_info = []
        self.layer_name_to_index = {}
        
        for i in range(self.layer_count):
            layer = self.dataset.GetLayer(i)
            layer_name = layer.GetName()

            self.layers.append(layer)
            self.layers_info.append({
                'name': layer_name,
                'features': len(layer),
                'fields': layer.GetLayerDefn().GetFieldCount(),
                'geomtypes': {ogr.GeometryTypeToName(k):v for k,v in layer.GetGeometryTypes().items()},
                'extent': layer.GetExtent(),
            })
            self.layer_name_to_index[layer_name] = i

    def print_info(self, layer_info):
        """
        prints layer information using the layer_info dict given.

        Parameters
        --------------------
        layer_info : dict
            dict containing 'name', 'features', 'fields', 'geomtypes', and 'extent' items
        """
        print("{} layer info:".format(layer_info['name']))
        print("features: {}".format(layer_info['features']))
        print("fields: {}".format(layer_info['fields']))
        for k,v in layer_info['geomtypes'].items():
            print("{} features of geomtype {}".format(v, k))
        print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(*layer_info['extent']))
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
            for layer_info in self.layers_info:
                self.print_info(layer_info)
        elif type(layer) == str:
            print_info(layers_info[layer_name_to_index[layer]])
        elif type(layer) == int:
            print_info(layers_info[layer])
        else:
              TypeError("layer parameter cannot be of type {}".format(type(layer)))

'''
TODO:
 - add a plotting function using matplotlib (specifically plt.Polygon and add patch)
    stackoverflow.com/questions/30447790
'''
