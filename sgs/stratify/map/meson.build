if WIN
  #create breaks module
  py.extension_module(
    'map_stratifications',
    'map_stratifications.cpp',
    link_with: [
      raster_lib,
    ],
    link_args: win_link_args,
    dependencies: [lib_deps],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    install: true,
  )
else
  #create srs module
  py.extension_module(
    'map_stratifications',
    'map_stratifications.cpp',
    dependencies: [lib_deps],
    include_directories: [
      inc_dir,
      utils_inc_dir,
    ],
    objects: [
      #from sgs/utils
      raster_obj,
    ],
    install: true,
  )
endif
