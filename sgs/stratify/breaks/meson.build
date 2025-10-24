#create breaks module
py.extension_module(
  'breaks',
  'breaks.cpp',
  link_args: win_link_args,
  dependencies: [
    lib_deps, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir,
  ],
  install: true,
)
