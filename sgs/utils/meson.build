#create raster module
py.extension_module(
  'raster',
  ['raster_pybind.cpp'],
  link_args: link_args, 
  dependencies: [
    gdal_dep,
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir
  ],
  install: true,
  subdir: 'sgspy',
)

#create vector module
py.extension_module(
  'vector',
  ['vector_pybind.cpp'],
  link_args: link_args,
  dependencies: [
    gdal_dep, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir
  ],
  install: true,
  subdir: 'sgspy',
)
