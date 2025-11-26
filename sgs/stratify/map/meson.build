#create map module
py.extension_module(
  'map_stratifications',
  'map_stratifications.cpp',
  link_args: link_args,
  dependencies: [
    gdal_dep, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir,
  ],
  install: true,
)
