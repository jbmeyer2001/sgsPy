#create quantiles module
py.extension_module(
  'quantiles',
  'quantiles.cpp',
  link_args: link_args,
  dependencies: [
    gdal_dep,
    mkl_dep,
    tbb_dep,
    pybind11_dep
  ],
  include_directories: [
    pip_inc_dir,
    inc_dir,
    utils_inc_dir,
  ],
  cpp_args: [
    '-DMKL_ILP64',
    '-m64'
  ],
  install: true,
)
