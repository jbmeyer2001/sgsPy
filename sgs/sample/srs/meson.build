#create srs module
py.extension_module(
  'srs',
  'srs.cpp',
  link_args: link_args,
  dependencies: [
    gdal_dep, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir,
    xoshiro_inc_dir,
  ],
  install: true,
)
