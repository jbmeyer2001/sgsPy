project(
  'sgs',
  'cpp',
  default_options: [
    'cpp_std=c++20'
  ],
)

#get build operating system
WIN = host_machine.system() == 'windows' or host_machine.system() == 'cygwin'

#get python installation
py = import('python').find_installation(pure: false)

#get compiler
cxx = meson.get_compiler('cpp')

#this is some very hacky build code... 
#
#The problem is that I have a build dependency (daal). Trying to compile daal
#from source opens up a massive massive can of worms.
#
#The daal python module, for some reason, doesn't have pkg-config or cmake information
#anywhere, so it isn't discoverable as a dependency by meson using the dependency()
#function (like pybind11 is).
#
#Thus, the only way to access those headers and libraries is to find the path of the 
#isolated build environment created by pip. This directory isn't given to meson directly --
#rather, the directory is prepended to the PATH environment variable before calling any
#meson build commands.
#
#meson doesn't have direct access to this variable (probably for good reason...) so the way
#to actually get access to it is by running a quick python script which prints the PATH
#variable, and parsing the printed string in meson...
if WIN
  path_entries = run_command(py, '-c', 'import os; print(os.environ["PATH"])', check: true).stdout().strip().split(';')
  foreach entry : path_entries
    if entry.contains('Temp') and entry.contains('pip-build-env') and entry.contains('overlay') and entry.contains('Scripts')
      pip_build_dir = entry.split('Scripts').get(0)
      pip_inc_dir = pip_build_dir + join_paths('Library', 'include')
      pip_lib_dir = pip_build_dir + join_paths('Library', 'lib')
      pip_bin_dir = pip_build_dir + join_paths('Library', 'bin')
    endif
  endforeach
else
  pip_build_dir = run_command(py, '-c', 'import os; print(os.environ["PATH"])', check: true).stdout().strip().split(':').get(0).split('bin').get(0)

  pip_inc_dir = pip_build_dir + 'include'
  pip_lib_dir = pip_build_dir + 'lib'
  pip_bin_dir = pip_build_dir + 'bin'
endif

#define directories
if WIN
  triplet = 'x64-windows'
else #linux
  triplet = 'x64-linux'
endif
vcpkg_inc_dir = include_directories('sgs/extern/vcpkg/installed/' + triplet + '/include')
xoshiro_inc_dir = include_directories('sgs/extern/xoshiro/include')
lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/' + triplet + '/lib'
bin_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/' + triplet + '/bin'
inc_dir = include_directories('sgs')

#get gdal dependency from vcpkg pkg-config
gdal_dep = dependency('gdal', method: 'pkg-config')

#get pybind11 dep from pyproject.toml build dependencies
pybind11_dep = dependency('pybind11')

#specify libs depending on os
if WIN
  onedal_libs = ['onedal_dll', 'onedal_core_dll']
  tbb_libs = ['tbb12']
  mkl_libs = ['mkl_intel_ilp64_dll', 'mkl_tbb_thread_dll', 'mkl_core_dll']
else #linux
  onedal_libs = ['onedal', 'onedal_core', 'onedal_parameters', 'onedal_thread']
  tbb_libs = ['tbb']
  mkl_libs = ['mkl_intel_ilp64', 'mkl_tbb_thread', 'mkl_core']

  #some libs are missing their non-versioned simlinks. Create them for meson to find.
  simlinks = [
    ['onedal', '.3'], 
    ['onedal_core', '.3'],
    ['onedal_parameters', '.3'], 
    ['onedal_thread', '.3'], 
    ['mkl_intel_ilp64', '.2'],
    ['mkl_tbb_thread', '.2'],
    ['mkl_core', '.2'],
  ]
  foreach simlink : simlinks
    #manually create simlink
    lib = simlink.get(0)
    version = simlink.get(1)
    
    simlink_src = pip_lib_dir + '/lib' + lib + '.so'
    simlink_dst = pip_lib_dir + '/lib' + lib + '.so' + version #this is the file which actually exists
    
    run_command('ln', '-s', simlink_dst, simlink_src, check: true)
  endforeach
endif

#add libs for oneDAL build dependency
oneDAL_dep = []
foreach lib : onedal_libs
  oneDAL_dep += cxx.find_library(lib, dirs: [pip_lib_dir])
endforeach
 
#add libs for tbb build dependency
tbb_dep = []
foreach lib : tbb_libs
  tbb_dep += cxx.find_library(lib, dirs: [pip_lib_dir])
endforeach

#add libs for mkl build dependency
mkl_dep = []
foreach lib : mkl_libs
  mkl_dep += cxx.find_library(lib, dirs: [pip_lib_dir])
endforeach

#add link args
link_args = []
if WIN
  link_args = [
    '-lws2_32',   #ws2_32.lib from Windows SDK 
    '-lcrypt32',  #crypt32.lib from Windows SDK
    '-lsecur32',  #secur32.lib from Windows SDK
    '-lwldap32',  #wldap32.lib from Windows SDK
    '-lwbemuuid'  #wbemuuid.lib from Windows SDK
  ]
else
  link_args = [
    '-Wl,-rpath,$ORIGIN/../../../sgs',
    '-Wl,-rpath,$ORIGIN/../../',
  ]
endif
add_global_link_arguments('-L' + lib_dir, language: 'cpp')

if WIN
  required_files = [
    #vcpkg dlls
    join_paths(bin_dir, 'aec.dll'),
    join_paths(bin_dir, 'bz2.dll'),
    join_paths(bin_dir, 'charset-1.dll'),
    join_paths(bin_dir, 'ffi-8.dll'),
    join_paths(bin_dir, 'freexl-1.dll'),
    join_paths(bin_dir, 'gdal.dll'),
    join_paths(bin_dir, 'geos.dll'),
    join_paths(bin_dir, 'geos_c.dll'),
    join_paths(bin_dir, 'geotiff.dll'),
    join_paths(bin_dir, 'gif.dll'),
    join_paths(bin_dir, 'hdf5.dll'),
    join_paths(bin_dir, 'hdf5_cpp.dll'),
    join_paths(bin_dir, 'hdf5_hl.dll'),
    join_paths(bin_dir, 'hdf5_hl_cpp.dll'),
    join_paths(bin_dir, 'iconv-2.dll'),
    join_paths(bin_dir, 'jpeg62.dll'),
    join_paths(bin_dir, 'json-c.dll'),
    join_paths(bin_dir, 'legacy.dll'),
    join_paths(bin_dir, 'Lerc.dll'),
    join_paths(bin_dir, 'libcrypto-3-x64.dll'),
    join_paths(bin_dir, 'libcurl.dll'),
    join_paths(bin_dir, 'libecpg.dll'),
    join_paths(bin_dir, 'libecpg_compat.dll'),
    join_paths(bin_dir, 'libexpat.dll'),
    join_paths(bin_dir, 'liblzma.dll'),
    join_paths(bin_dir, 'libpgtypes.dll'),
    join_paths(bin_dir, 'libpng16.dll'),
    join_paths(bin_dir, 'libpq.dll'),
    join_paths(bin_dir, 'libsharpyuv.dll'),
    join_paths(bin_dir, 'libssl-3-x64.dll'),
    join_paths(bin_dir, 'libwebp.dll'),
    join_paths(bin_dir, 'libwebpdecoder.dll'),
    join_paths(bin_dir, 'libwebpdemux.dll'),
    join_paths(bin_dir, 'libwebpmux.dll'),
    join_paths(bin_dir, 'libxml2.dll'),
    join_paths(bin_dir, 'lz4.dll'),
    join_paths(bin_dir, 'minizip.dll'),
    join_paths(bin_dir, 'netcdf.dll'),
    join_paths(bin_dir, 'openjp2.dll'),
    join_paths(bin_dir, 'pcre2-16.dll'),
    join_paths(bin_dir, 'pcre2-32.dll'),
    join_paths(bin_dir, 'pcre2-8.dll'),
    join_paths(bin_dir, 'pcre2-posix.dll'),
    join_paths(bin_dir, 'pkgconf-7.dll'),
    join_paths(bin_dir, 'proj_9.dll'),
    join_paths(bin_dir, 'qhull_r.dll'),
    join_paths(bin_dir, 'spatialite.dll'),
    join_paths(bin_dir, 'sqlite3.dll'),
    join_paths(bin_dir, 'szip.dll'),
    join_paths(bin_dir, 'tiff.dll'),
    join_paths(bin_dir, 'turbojpeg.dll'),
    join_paths(bin_dir, 'uriparser.dll'),
    join_paths(bin_dir, 'zlib1.dll'),
    join_paths(bin_dir, 'zstd.dll'),
  
    #onedal dlls
    join_paths(pip_bin_dir, 'onedal.3.dll'),
    join_paths(pip_bin_dir, 'onedal_core.3.dll'),
    join_paths(pip_bin_dir, 'onedal_thread.3.dll'),
    join_paths(pip_bin_dir, 'tbb12.dll'),
    join_paths(pip_bin_dir, 'tbbmalloc.dll'),
  
    #mkl dlls
    join_paths(pip_bin_dir, 'mkl_core.2.dll'),
    join_paths(pip_bin_dir, 'mkl_tbb_thread.2.dll'),
    join_paths(pip_bin_dir, 'mkl_blacs_intelmpi_ilp64.2.dll'),
    join_paths(pip_bin_dir, 'mkl_vml_def.2.dll'),
    join_paths(pip_bin_dir, 'mkl_vml_avx2.2.dll'),
    join_paths(pip_bin_dir, 'mkl_vml_avx512.2.dll'),
    join_paths(pip_bin_dir, 'mkl_vml_cmpt.2.dll'),
    join_paths(pip_bin_dir, 'mkl_vml_mc3.2.dll'),

    #proj
    meson.current_source_dir() + '/sgs/extern/vcpkg/installed/' + triplet + '/share/proj/proj.db'
  ]
else #linux
  required_files = []
  #only vendor the files which need their rpath to be specifically edited
  vendored_libs = [
    #onedal libs
    join_paths(pip_lib_dir, 'libonedal.so.3'),
    #join_paths(pip_lib_dir, 'libonedal_core.so.3'),
    #join_paths(pip_lib_dir, 'libonedal_parameters.so.3'),
    #join_paths(pip_lib_dir, 'libonedal_thread.so.3'),
    #join_paths(pip_lib_dir, 'libtbb.so'),
    #join_paths(pip_lib_dir, 'libtbbmalloc.so'),
    #join_paths(pip_lib_dir, 'libtbb.so.12'),
    #join_paths(pip_lib_dir, 'libtbbmalloc.so.2'),
    #join_paths(pip_lib_dir, 'libtbb.so.12.17'),
    #join_paths(pip_lib_dir, 'libtbbmalloc.so.2.17'),
  
    #mkl libs
    #join_paths(pip_lib_dir, 'libmkl_intel_ilp64.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_tbb_thread.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_core.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_vml_def.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_vml_avx2.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_vml_avx512.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_vml_cmpt.so.2'),
    #join_paths(pip_lib_dir, 'libmkl_vml_mc3.so.2'),

    #proj
    #meson.current_source_dir() + '/sgs/extern/vcpkg/installed/' + triplet + '/share/proj/proj.db'
  ]
endif

foreach file : required_files
  install_data(file)
endforeach

foreach lib_path : vendored_libs
  run_command('patchelf', '--set-rpath', '$ORIGIN/../lib', lib_path)
  install_data(lib_path)
endforeach

py.extension_module(
  '_sgs',
  'sgs/pybind.cpp',
  link_args: link_args,
  dependencies: [
    gdal_dep,
    oneDAL_dep,
    mkl_dep,
    tbb_dep,
    pybind11_dep
  ],
  include_directories: [
    pip_inc_dir,
    vcpkg_inc_dir,
    xoshiro_inc_dir,
    inc_dir
  ],
  cpp_args: [
    '-DMKL_ILP64',
    '-m64'
  ],
  install: true,
)

