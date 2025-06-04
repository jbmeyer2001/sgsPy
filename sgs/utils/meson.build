#raster_lib = shared_library(
#'raster_lib', 
#sources: ['raster.cpp', 'raster.h'],
#dependencies: [lib_deps],
#include_directories: inc_dir,
#)

py.extension_module(
  'raster',
  ['raster.cpp', 'raster.h', 'raster_pybind.cpp'],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  install: true,
)

py.extension_module(
  'vector',
  ['vector.cpp', 'vector.h', 'vector_pybind.cpp'],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  install: true,
)
