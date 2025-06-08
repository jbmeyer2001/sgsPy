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

#both windows and linux builds need to be linked slightly differently.
#windows requires the use of the link_with argument, whereas Linux
#requires extracted objects to be added.
#
#Windows additionally requires linker arguments to search for .lib files
#from the Windows SDK.
if WIN
  #create raster module
  py.extension_module(
    'raster',
    ['raster_pybind.cpp'],
    link_with: [raster_lib],
    link_args: win_link_args, 
    dependencies: [lib_deps],
    include_directories: inc_dir,
    install: true,
  )

  #create vector module
  py.extension_module(
    'vector',
    ['vector_pybind.cpp'],
    link_with: [vector_lib],
    link_args: win_link_args,
    dependencies: [lib_deps],
    include_directories: inc_dir,
    install: true,
  )
else #linux
  #first, extract objects to be linked. These will be used elsewhere in the project as well
  raster_obj = raster_lib.extract_objects('raster.cpp')
  vector_obj = vector_lib.extract_objects('vector.cpp')

  #create raster module
  #py.extension_module(
    #'raster',
    #['raster_pybind.cpp'],
    #dependencies: [lib_deps],
    #include_directories: inc_dir,
    #objects: [raster_obj],
    #install: true,
    #)

  #create vector module
  #py.extension_module(
    #'vector',
    #['vector_pybind.cpp'],
    #dependencies: [lib_deps],
    #include_directories: inc_dir,
    #objects: [vector_obj],
    #install: true,
    #  )
endif
