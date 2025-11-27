#create systematic module
py.extension_module(
  'systematic',
  'systematic.cpp',
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
