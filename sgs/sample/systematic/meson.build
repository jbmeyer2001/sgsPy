#if WIN
  #create systematic module
  py.extension_module(
    'systematic',
    'systematic.cpp',
    dependencies: [
      lib_deps,
      pybind11_dep
    ],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    install: true,
    link_args: win_link_args,
  )

