if WIN
  #create breaks module
  py.extension_module(
    'pca', 'pca.cpp',
    #input_files,
    link_args: win_link_args,
    dependencies: [lib_deps, tbb_lib_deps, oneDAL_lib_deps, pybind11_dep],
    include_directories: [
      inc_dir,
      utils_inc_dir,
      pip_inc_dir,
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
