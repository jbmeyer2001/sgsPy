if WIN
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
    install: true,
  )
else
  #create quantiles module
  py.extension_module(
    'quantiles',
    'quantiles.cpp',
    dependencies: [lib_deps, mkl_lib_deps, tbb_lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    install: true,
  )
endif
