#both windows and linux builds need to be linked slightly differently.
#windows requires the use of the link_with argument, whereas Linux
#requires extracted objects to be added.
#
#Windows additionally requires linker arguments to search for .lib files
#from the Windows SDK.

if WIN
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
    install: true,
    link_args: win_link_args,
  )
else #linux
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
      #from sgs/utils
      raster_obj,
      vector_obj,

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
    ],
    install: true,
)
endif

