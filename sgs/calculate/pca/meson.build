#create pca module
py.extension_module(
  'pca', 'pca.cpp',
  link_args: link_args,
  #input_files,
  dependencies: [
    gdal_dep, 
    oneDAL_lib_deps, 
    pybind11_dep
  ],
  include_directories: [
    inc_dir,
    utils_inc_dir,
    pip_inc_dir,
  ],
  install: true,
)
