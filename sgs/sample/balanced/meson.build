#create balanced module
py.extension_module(
  'balanced',
  ['balanced.cpp'],
  dependencies: [lib_deps],
  include_directories: [inc_dir, utils_inc_dir],
  objects: [raster_obj, vector_obj],
  install: true,
)
