#create balanced module
py.extension_module(
  'balanced',
  ['balanced.cpp'],
  dependencies: [lib_deps],
  include_directories: [
    inc_dir, 
    utils_inc_dir, 
    balanced_sampling_inc_dir
  ],
  objects: [
    raster_obj,                 #from sgs/utils
    vector_obj,                 #from sgs/utils
    cube_class_obj,             #from sgs/extern/BalancedSampling
    cube_stratified_class_obj,  #form sgs/extern/BalancedSampling
    index_list_class_obj,       #from sgs/extern/BalancedSampling
    kdtree_class_obj,           #from sgs/extern/BalancedSampling
    lpm_class_obj               #from sgs/extern/BalancedSampling
  ],
  install: true,
)
