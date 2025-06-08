#build GDALRasterWrapper as library
raster_lib = static_library(
  'raster_lib', 'raster.cpp', 'raster.h',
  dependencies: [lib_deps, py.dependency()],
  include_directories: [inc_dir],
)

#build GDALVectorWrapper as library
vector_lib = static_library(
  'vector_lib', 'vector.cpp', 'vector.h',
  dependencies: [lib_deps, py.dependency()],
  include_directories: [inc_dir],
)

#acquire objects if necessary
if not WIN
  raster_obj = raster_lib.extract_objects('raster.cpp')
  vector_obj = vector_lib.extract_objects('vector.cpp')
  raster_module_objects = [raster_obj]
  vector_module_objects = [vector_obj]
else
  raster_module_objects = []
  vector_module_objects = []
endif

#create raster module
py.extension_module(
  'raster',
  ['raster_pybind.cpp'],
  link_with: [raster_lib],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  objects: raster_module_objects, #only defined if linux
  install: true,
  link_args: link_arguments,    #only defined if windows
)

#create vector module
py.extension_module(
  'vector',
  ['vector_pybind.cpp'],
  link_with: [vector_lib],
  dependencies: [lib_deps],
  include_directories: inc_dir,
  objects: vector_module_objects, #only defined if linux
  install: true,
  link_args: link_arguments,      #only defined if windows
)
