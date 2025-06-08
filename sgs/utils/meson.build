#build GDALRasterWrapper as library

#path_to_install = py.get_install_dir() + 'sgs/utils'
#message(path_to_install)

#raster_lib = shared_library(
#  'raster', 'raster.cpp', 'raster.h',
#  dependencies: [lib_deps, py.dependency()],
#  include_directories: [inc_dir],
#  link_args: [ '/IMPLIB:' + join_paths(meson.current_build_dir(), 'sgs/utils/raster.lib') ],
#  install: true,
#)

#build GDALVectorWrapper as library
#vector_lib = shared_library(
#  'vector_lib', 'vector.cpp', 'vector.h',
#  dependencies: [lib_deps, py.dependency()],
#  include_directories: [inc_dir],
#  install: true
#)

#acquire objects to be used in python modules
#raster_obj = raster_lib.extract_objects('raster.cpp')
#vector_obj = vector_lib.extract_objects('vector.cpp')

#create raster module
py.extension_module(
  'raster',
  ['raster.h', 'raster.cpp', 'raster_pybind.cpp'],
  #  link_with: [raster_lib],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  install: true,
)

#create vector module
#py.extension_module(
  #'vector',
  #['vector_pybind.cpp'],
  #dependencies: [lib_deps],
  #include_directories: inc_dir,
  #  objects: [vector_obj],
#  install: true,
  #)
