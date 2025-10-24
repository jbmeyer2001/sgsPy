if WIN
  #create breaks module
  py.extension_module(
    'pca', 'pca.cpp',
    #input_files,
    link_with: [
      raster_lib,
    ],
    link_args: win_link_args,
    dependencies: [lib_deps],
    include_directories: [
      inc_dir,
      utils_inc_dir,
      oneapi_inc_dir,
      oneapi_daal_inc_dir,
      oneapi_daal_include_inc_dir
    ],
    cpp_args: [
      '-D__ONEDAL_IDE_MODE__', #read why I do this here: 
                              #https://github.com/uxlfoundation/oneDAL/blob/340926e2284d2c1cc87238538fcba2b265aa4980/cpp/oneapi/dal/backend/dispatcher_cpu.hpp#L22
    ],
    install: true,
  )
else
  #create srs module
  py.extension_module(
    'pca', 'pca.cpp',
    #input_files,
    dependencies: [lib_deps, tbb_lib_deps, oneDAL_lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
      pip_inc_dir,
    ],
    install: true,
  )
endif
