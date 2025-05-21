from osgeo import gdal, ogr
import matplotlib.pyplot as ptl

class SpatialVector:
    def __init__(self, image):
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
        print("{} layer info:".format(layer_info['name']))
        print("features: {}".format(layer_info['features']))
        print("fields: {}".format(layer_info['fields']))
        for k,v in layer_info['geomtypes'].items():
            print("{} features of geomtype {}".format(v, k))
        print("bounds (xmin, xmax, ymin, ymax): ({}, {}, {}, {})".format(*layer_info['extent']))
        print()

    def info(self, layer=None):
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
