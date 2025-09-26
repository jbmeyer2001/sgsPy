#build GDALRasterWrapper as library
#raster_lib = library(
#  'raster_lib', 'raster.cpp', 'raster.h',
#  dependencies: [lib_deps, py.dependency()],
#  include_directories: [inc_dir],
#)

#raster_lib_dep = declare_dependency(
#  link_with: raster_lib,
#  include_directories: [inc_dir, utils_inc_dir]
#)

#build GDALVectorWrapper as library
#vector_lib = library(
#  'vector_lib', 'vector.cpp', 'vector.h',
#  dependencies: [lib_deps, py.dependency()],
#  include_directories: [inc_dir],
#)

#vector_lib_dep = declare_dependency(
#  link_with: vector_lib,
#  include_directories: [inc_dir, utils_inc_dir]
#)

#both windows and linux builds need to be linked slightly differently.
#windows requires the use of the link_with argument, whereas Linux
#requires extracted objects to be added.
#
#Windows additionally requires linker arguments to search for .lib files
#from the Windows SDK.
if WIN
  #create access lib
  #access_lib = library(
  #  'access_lib', 'access.cpp', 'access.h', 'vector.h',
  #  link_with: [raster_lib],
  #  link_args: win_link_args,
  #  dependencies: [lib_deps, py.dependency()],
  #  include_directories: inc_dir,
  #)

  #access_lib_dep = declare_dependency(
  #  link_with: access_lib,
  #  include_directories: [inc_dir, utils_inc_dir]
  #)

  #create raster module
  py.extension_module(
    'raster',
    ['raster_pybind.cpp', 'raster.h'],
    link_args: win_link_args, 
    dependencies: [lib_deps],
    include_directories: [inc_dir],
    install: true,
  )

  #create vector module
  py.extension_module(
    'vector',
    ['vector_pybind.cpp', 'vector.h'],
    link_args: win_link_args,
    dependencies: [lib_deps],
    include_directories: [inc_dir],
    install: true,
  )

else #linux
  #first, extract objects to be linked. These will be used elsewhere in the project as well
  raster_obj = raster_lib.extract_objects('raster.cpp')
  vector_obj = vector_lib.extract_objects('vector.cpp')

  #create access library and extract object, using raster and vector objects
  access_lib = library(
    'access_lib', 'access.cpp', 'access.h',
    dependencies: [lib_deps, py.dependency()],
    include_directories: inc_dir,
    objects: [raster_obj, vector_obj],
  )
  access_obj = access_lib.extract_objects('access.cpp')

  #create raster module
  py.extension_module(
    'raster',
    ['raster_pybind.cpp'],
    dependencies: [lib_deps],
    include_directories: inc_dir,
    objects: [raster_obj],
    install: true,
  )

  #create vector module
  py.extension_module(
    'vector',
    ['vector_pybind.cpp'],
    dependencies: [lib_deps, py.dependency()],
    include_directories: inc_dir,
    objects: [vector_obj],
    install: true,
  )
endif
