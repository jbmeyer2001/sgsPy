daal_pca_dir = root_dir + '/sgs/extern/oneDAL/cpp/daal/src/algorithms/pca/'
daal_data_management_dir = root_dir + '/sgs/extern/oneDAL/cpp/daal/src/data_management/'
input_files = [
  'pca.cpp',
  daal_data_management_dir + 'homogen_numeric_table.cpp',
  daal_pca_dir + 'pca_dense_correlation_distr_step2_kernel.h',       	  
  daal_pca_dir + 'pca_dense_svd_online_kernel.h',
  daal_pca_dir + 'pca_baseparameter.cpp',                                 
  daal_pca_dir + 'pca_dense_correlation_online_container.h',         
  daal_pca_dir + 'pca_distributedinput_correlation.cpp',
  daal_pca_dir + 'pca_baseparameter_fpt.cpp',                        
  daal_pca_dir + 'pca_dense_correlation_online_fpt_cpu.cpp',         
  daal_pca_dir + 'pca_distributedinput_svd.cpp',
  daal_pca_dir + 'pca_batchparameter_correlation_fpt.cpp',              
  daal_pca_dir + 'pca_dense_correlation_online_fpt_dispatcher.cpp',  
  daal_pca_dir + 'pca_distributedparameter_correlation.h',
  daal_pca_dir + 'pca_batchparameter_svd_fpt.cpp',                              
  daal_pca_dir + 'pca_distributedparameter_correlation_fpt.cpp',
  daal_pca_dir + 'pca_dense_base.h',                                    
  daal_pca_dir + 'pca_dense_correlation_online_kernel.h',            
  daal_pca_dir + 'pca_input.cpp',
  daal_pca_dir + 'pca_dense_base_fpt_cpu.cpp',                          
  daal_pca_dir + 'pca_dense_svd_base.h',                             
  daal_pca_dir + 'pca_onlineparameter_correlation.h',
  daal_pca_dir + 'pca_dense_base_iface.h',                            
  daal_pca_dir + 'pca_dense_svd_batch_container.h',                  
  daal_pca_dir + 'pca_onlineparameter_correlation_fpt.cpp',                            
  daal_pca_dir + 'pca_dense_svd_batch_fpt_cpu.cpp',                  
  daal_pca_dir + 'pca_onlineparameter_svd.h',
  daal_pca_dir + 'pca_dense_correlation_base.h',                     
  daal_pca_dir + 'pca_dense_svd_batch_fpt_dispatcher.cpp',           
  daal_pca_dir + 'pca_onlineparameter_svd_fpt.cpp',
  daal_pca_dir + 'pca_dense_correlation_base_fpt_cpu.cpp',                               
  daal_pca_dir + 'pca_partialresult_correlation.cpp',
  daal_pca_dir + 'pca_dense_correlation_base_iface.h',                
  daal_pca_dir + 'pca_dense_svd_batch_kernel.h',                     
  daal_pca_dir + 'pca_partialresult_correlation.h',             
  daal_pca_dir + 'pca_dense_svd_batch_kernel_instance.h',            
  daal_pca_dir + 'pca_partialresult_correlation_fpt.cpp',
  daal_pca_dir + 'pca_dense_correlation_batch_container.h',          
  daal_pca_dir + 'pca_dense_svd_container.h',                        
  daal_pca_dir + 'pca_partialresult_svd.cpp',
  daal_pca_dir + 'pca_dense_correlation_batch_fpt_cpu.cpp',           
  daal_pca_dir + 'pca_dense_svd_distr_step2_container.h',            
  daal_pca_dir + 'pca_partialresult_svd.h',
  daal_pca_dir + 'pca_dense_correlation_batch_fpt_dispatcher.cpp',     
  daal_pca_dir + 'pca_dense_svd_distr_step2_fpt_cpu.cpp',            
  daal_pca_dir + 'pca_partialresult_svd_fpt.cpp',             
  daal_pca_dir + 'pca_dense_svd_distr_step2_fpt_dispatcher.cpp',     
  daal_pca_dir + 'pca_result.cpp',
  daal_pca_dir + 'pca_dense_correlation_batch_kernel.h',                           
  daal_pca_dir + 'pca_result_fpt.cpp',
  daal_pca_dir + 'pca_dense_correlation_batch_kernel_instance.h',      
  daal_pca_dir + 'pca_dense_svd_distr_step2_kernel.h',               
  daal_pca_dir + 'pca_result_impl.cpp',
  daal_pca_dir + 'pca_dense_correlation_distr_step2_container.h',       
  daal_pca_dir + 'pca_dense_svd_online_container.h',                 
  daal_pca_dir + 'pca_result_impl.h',
  daal_pca_dir + 'pca_dense_correlation_distr_step2_fpt_cpu.cpp',    
  daal_pca_dir + 'pca_dense_svd_online_fpt_cpu.cpp',                 
  daal_pca_dir + 'pca_result_impl_fpt.cpp',
  daal_pca_dir + 'pca_dense_correlation_distr_step2_fpt_dispatcher.cpp',  
  daal_pca_dir + 'pca_dense_svd_online_fpt_dispatcher.cpp'
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
      daal_inc_dir
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
      daal_inc_dir
    ],
    objects: [
      #from sgs/utils
      raster_obj,
    ],
    install: true,
  )
endif
