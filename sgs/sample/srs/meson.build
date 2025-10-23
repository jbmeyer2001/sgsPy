if WIN
  #create srs module
  py.extension_module(
    'srs',
    'srs.cpp',
    link_with: [
      raster_lib,
      vector_lib,
      access_lib,
    ],
    dependencies: [lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    install: true,
    link_args: win_link_args,
  )
else #linux
  #create srs module
  py.extension_module(
    'srs',
    'srs.cpp',
    dependencies: [lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    objects: [
      #from sgs/utils
      raster_obj,
      vector_obj,
      access_obj,
    ],
    install: true,
  )
endif
