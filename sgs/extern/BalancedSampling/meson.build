#build balanced sampling package as library
balanced_sampling_lib = static_library(
  'balanced_sampling_lib',
  'CubeClass.cc',
  'CubeClass.h',
  'CubeStratifiedClass.cc',
  'CubeStratifiedClass.h',
  'CubeVectorInNullSpace.cc',
  'CubeVectorInNullSpace.h',
  'IndexListClass.cc',
  'IndexListClass.h',
  'KDNodeClass.cc',
  'KDNodeClass.h',
  'KDStoreClass.cc',
  'KDStoreClass.h',
  'KDTreeClass.cc',
  'KDTreeClass.h',
  'LpmClass.cc',
  'LpmClass.h',
  'ReducedRowEchelonForm.cc',
  'ReducedRowEchelonForm.h',
  'uniform.h',
  'utils.h',
  'utils-matrix.h',
)

#the classes which are included directly in balanced.cpp
cube_class_obj = balanced_sampling_lib.extract_objects('CubeClass.cc')
cube_stratified_class_obj = balanced_sampling_lib.extract_objects('CubeStratifiedClass.cc')
index_list_class_obj = balanced_sampling_lib.extract_objects('IndexListClass.cc')
kdtree_class_obj = balanced_sampling_lib.extract_objects('KDTreeClass.cc')
lpm_class_obj = balanced_sampling_lib.extract_objects('LpmClass.cc')
