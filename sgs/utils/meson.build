raster_lib = static_library(
  'raster_lib', 'raster.cpp', 'raster.h',
  dependencies: [lib_deps, py.dependency()],
  include_directories: [inc_dir],
)

raster_obj = raster_lib.extract_objects('raster.cpp')

py.extension_module(
  'raster',
  ['raster_pybind.cpp'],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  objects: [raster_obj],
  install: true,
)

py.extension_module(
  'vector',
  ['vector.cpp', 'vector.h', 'vector_pybind.cpp'],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  install: true,
)
