#create raster module
py.extension_module(
  'raster',
  ['raster_pybind.cpp'],
  link_args: win_link_args, 
  dependencies: [
    lib_deps,
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir
  ],
  install: true,
)

#create vector module
py.extension_module(
  'vector',
  ['vector_pybind.cpp'],
  link_args: win_link_args,
  dependencies: [
    lib_deps, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir
  ],
  install: true,
)
