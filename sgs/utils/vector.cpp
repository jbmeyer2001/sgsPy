#include <gdal_priv.h>
#include <gdal.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

/**
 *
 */
class SpatialVector {
	private:
	typedef struct layerInfo {
		std::unique_ptr<OGRLayer> layer;
		std::string name;
		int features;
		int fields;
		std::unordered_map<std::string, int>;
		double extent[4];
	};
	std::unique_ptr<GDALDataset> p_dataset;
	int numLayers;
	std::vector<layerInfo> layers;
	std::unordered_map<std::string, int> layerNameMap;

	public:
	
	/**
	 *
	 */	
	SpatialVector(std::string filename) {
		//dataset
		this->p_dataset = std::unique_ptr<GDALDataset>(GDALDataset::FromHandle(GDALOpen(filename.c_str(), GA_ReadOnly)));

		//layers
		this->numLayers = this->p_dataset->GetLayerCount();
		for (int i = 0; i < this->numLayers; i++) {
			layerInfo newLayer;
			
			newLayer.layer = std::unique_ptr<OGRLayer>(this->p_dataset->GetLayer(i));
			newLayer.name = std::string(newLayer.layer->GetName());
			
			int featureCount	
		}
	}
};
