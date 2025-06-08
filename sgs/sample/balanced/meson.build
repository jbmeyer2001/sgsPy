#create balanced module
py.extension_module(
  'balanced',
  ['balanced.cpp'],
  link_with: [
    raster_lib,
    vector_lib,
    balanced_sampling_lib,
  ],
  dependencies: [lib_deps],
  include_directories: [
    inc_dir, 
    utils_inc_dir, 
    balanced_sampling_inc_dir
  ],
  install: true,
  link_args: ['-lws2_32', '-lcrypt32', '-lsecur32', '-lwldap32', '-lwbemuuid'],
)
