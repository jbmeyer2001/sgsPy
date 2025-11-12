if WIN
  #create quantiles module
  py.extension_module(
    'quantiles',
    'quantiles.cpp',
    link_args: [
      win_link_args, 
    ],
    dependencies: [
      lib_deps,
      pybind11_dep
    ],
    link_args: [
      win_link_args, 
    ],
    dependencies: [
      lib_deps
    ],
    include_directories: [
      inc_dir,
      utils_inc_dir,
      mkl_inc_dir,
    ],
    cpp_args: [
      '-DMKL_ILP64',
      '-m64',
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
