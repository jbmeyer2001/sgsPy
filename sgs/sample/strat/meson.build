#create strat module
py.extension_module(
  'strat',
  'strat.cpp',
  link_args: win_link_args,
  dependencies: [
    lib_deps, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir,
    xoshiro_inc_dir
  ],
  install: true,
)
