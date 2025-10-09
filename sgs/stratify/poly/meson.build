if WIN
  #create breaks module
  py.extension_module(
    'poly',
    'poly.cpp',
    link_with: [
      raster_lib,
      vector_lib,
    ],
    link_args: win_link_args,
    dependencies: [lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    install: true,
  )
else
  #create srs module
  py.extension_module(
    'poly',
    'poly.cpp',
    dependencies: [lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    objects: [
      #from sgs/utils
      raster_obj,
      vector_obj,
    ],
    install: true,
  )
endif
