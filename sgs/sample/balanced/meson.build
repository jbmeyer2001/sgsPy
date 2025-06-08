#acquire objects if necessary
if not WIN
  sample_balanced_module_objects = [
    #from sgs/utils
    raster_lib,
    vector_lib,

    #from sgs/extern/BalancedSampling
    cube_class_obj,
    cube_stratified_class_obj,
    cube_vector_in_null_space_obj,
    index_list_class_obj,
    kdnode_class_obj,
    kdstore_class_obj,
    kdtree_class_obj,
    lpm_class_obj,
    rref_obj,
  ]
else
  sample_balanced_module_objects = []
endif

#create balanced module
py.extension_module(
  'balanced',
  ['balanced.cpp'],
  link_with: [
    raster_lib,
    vector_lib,
    balanced_sampling_lib,
  ],
  dependencies: [lib_deps],
  include_directories: [
    inc_dir, 
    utils_inc_dir, 
    balanced_sampling_inc_dir
  ],
  objects: sample_balanced_module_objects,  #only defined if linux
  install: true,
  link_args: link_arguments,                #only defined if windows
)
