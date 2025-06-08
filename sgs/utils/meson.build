#build GDALRasterWrapper as library
raster_lib = static_library(
  'raster_lib', 'raster.cpp', 'raster.h',
  dependencies: [lib_deps, py.dependency()],
  include_directories: [inc_dir],
)

#build GDALVectorWrapper as library
vector_lib = static_library(
  'vector_lib', 'vector.cpp', 'vector.h',
  dependencies: [lib_deps, py.dependency()],
  include_directories: [inc_dir],
)

#acquire objects to be used in python modules
#raster_obj = raster_lib.extract_objects('raster.cpp')
#vector_obj = vector_lib.extract_objects('vector.cpp')

#create raster module
py.extension_module(
  'raster',
  ['raster_pybind.cpp'],
  link_with: [raster_lib],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  #objects: [raster_obj],
  install: true,
  link_args: ['-lws2_32', '-lcrypt32', '-lsecur32', '-lwldap32', '-lwbemuuid'],
)

#create vector module
py.extension_module(
  'vector',
  ['vector_pybind.cpp'],
  link_with: [vector_lib],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  install: true,
  link_args: ['-lws2_32', '-lcrypt32', '-lsecur32', '-lwldap32', '-lwbemuuid'],
)
