oneapi_dir = root_dir + '/sgs/extern/oneDAL/cpp/oneapi/dal/'
input_files = [
  'pca.cpp',
  oneapi_dir + 'exceptions.cpp',
  oneapi_dir + 'algo/pca/common.cpp',
  oneapi_dir + 'algo/pca/train_types.cpp',
  oneapi_dir + 'backend/memory.cpp',
  oneapi_dir + 'backend/communicator.cpp',
  oneapi_dir + 'detail/array_impl.cpp',
  oneapi_dir + 'detail/chunked_array_impl.cpp',
  oneapi_dir + 'detail/common.cpp',
  oneapi_dir + 'detail/cpu.cpp',
  oneapi_dir + 'detail/error_messages.cpp',
  oneapi_dir + 'detail/global_context.cpp',
  oneapi_dir + 'detail/global_context_impl.cpp',
  oneapi_dir + 'detail/memory_impl_host.cpp',
  oneapi_dir + 'detail/policy.cpp',
  oneapi_dir + 'detail/serialization.cpp',
  oneapi_dir + 'detail/spmd_policy.cpp',
  oneapi_dir + 'spmd/communicator.cpp',
  oneapi_dir + 'table/common.cpp',
  oneapi_dir + 'table/csr.cpp',
  oneapi_dir + 'table/heterogen.cpp',
  oneapi_dir + 'table/homogen.cpp',
  oneapi_dir + 'table/backend/csr_kernels.cpp',
  oneapi_dir + 'table/backend/heterogen_kernels.cpp',
  oneapi_dir + 'table/backend/homogen_kernels.cpp',
  oneapi_dir + 'table/backend/convert/common.cpp',
  oneapi_dir + 'table/backend/convert/common_convert.cpp',
  oneapi_dir + 'table/backend/convert/copy_convert.cpp',
  oneapi_dir + 'table/detail/table_builder.cpp',
  oneapi_dir + 'table/detail/table_kinds.cpp',
  oneapi_dir + 'table/detail/table_utils.cpp',
]

if WIN
  #create breaks module
  py.extension_module(
    'pca',
    input_files,
    link_with: [
      raster_lib,
    ],
    link_args: win_link_args,
    dependencies: [lib_deps],
    include_directories: [
      inc_dir,
      utils_inc_dir,
      oneapi_inc_dir,
      oneapi_daal_inc_dir
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
    'pca',
    input_files,
    dependencies: [lib_deps],
    include_directories: [
      inc_dir,
      utils_inc_dir,
      oneapi_inc_dir,
      oneapi_daal_inc_dir
    ],
    objects: [
      #from sgs/utils
      raster_obj,
    ],
    cpp_args: [
      '-D__ONEDAL_IDE_MODE__', #read why I do this here: 
                              #https://github.com/uxlfoundation/oneDAL/blob/340926e2284d2c1cc87238538fcba2b265aa4980/cpp/oneapi/dal/backend/dispatcher_cpu.hpp#L22
    ],
    install: true,
  )
endif
