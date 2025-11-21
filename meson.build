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
inc_dir = include_directories('sgs/extern/vcpkg/installed/' + triplet + '/include')
lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/' + triplet + '/lib'
bin_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/' + triplet + '/bin'
utils_inc_dir = include_directories('sgs/utils')
xoshiro_inc_dir = include_directories('sgs/extern/xoshiro/include')

lib_deps = []
mkl_lib_deps = []
oneDAL_lib_deps = []
tbb_lib_deps = []

#add all libraries in the vcpkg lib directory to the project
#the following is hacky meson build code to add all the external libraries required for the project from the vcpkg install folder
if WIN
 else #linux
  #in linux list all of the libs rather than programatically finding them
  #because of multiple definition errors. 
  libs = [
    'aec',
    'crypto',
    'curl',
    'ecpg',
    'ecpg_compat',
    'expat',
    'freexl',
    'gdal',
    'geos',
    'geos_c',
    'geotiff',
    'gif',
    'hdf5',
    'hdf5_cpp',
    'hdf5_hl',
    'hdf5_hl_cpp',
    'jpeg',
    'json-c',
    'kmlbase',
    'kmlconvenience',
    'kmldom',
    'kmlengine',
    'kmlregionator',
    'kmlxsd',
    'Lerc',
    'lz4',
    'lzma',
    'minizip',
    'netcdf',
    'openjp2',
    'pcre2-8',
    'pcre2-16',
    'pcre2-32',
    'pcre2-posix',
    'pgcommon',
    'pgport',
    'pgtypes',
    'pkgconf',
    'png',
    'png16',
    'pq',
    'proj',
    'qhullcpp',
    'qhullstatic_r',
    'sharpyuv',
    'spatialite',
    'sqlite3',
    'ssl',
    'sz',
    'tiff',
    'turbojpeg',
    'uriparser',
    'webp',
    'webpdecoder',
    'webpdemux',
    'webpmux',
    'xml2',
    'z',
    'zstd',
    'mod_spatialite',
  ]

  foreach lib : libs
    lib_deps += cxx.find_library(lib, dirs : [lib_dir], static:true)
  endforeach

  mkl_libs = [
    'mkl_blacs_intelmpi_lp64',
    'mkl_blas95_lp64',
    'mkl_cdft_core',
    'mkl_core',
    'mkl_intel_lp64',
    'mkl_lapack95_lp64',
    'mkl_scalapack_lp64',
    'mkl_tbb_thread',
  ]

  foreach lib : mkl_libs
    mkl_lib_deps += cxx.find_library(lib, dirs : [lib_dir, intel_lib_dir], static:true)
  endforeach

  tbb_libs = [
    'tbb',
    'tbbmalloc'
  ]

  foreach lib : tbb_libs
    tbb_lib_deps += cxx.find_library(lib, dirs : [lib_dir, intel_lib_dir], static:true)
  endforeach

  oneDAL_libs = [
    'onedal',
    'onedal_core',
    'onedal_parameters',
    'onedal_thread'
  ]
 
  foreach lib : oneDAL_libs
    oneDAL_lib_deps += cxx.find_library(lib, dirs: [pip_lib_dir])
  endforeach


  oneDAL_libs = [
    'onedal',
    'onedal_core',
    'onedal_parameters',
    'onedal_thread'
  ]

  foreach lib : oneDAL_libs
    oneDAL_lib_deps += cxx.find_library(lib, dirs: [pip_lib_dir])
  endforeach
endif

#get gdal, tbb, mkl dependencies from vcpkg pkg-config
gdal_dep = dependency('gdal', method: 'pkg-config')

#get pybind11 dep from pyproject.toml build dependencies
pybind11_dep = dependency('pybind11')

#get compiler to add libs to
cxx = meson.get_compiler('cpp')

#add libs from oneDAL build dependency
oneDAL_dep = []
foreach lib : ['onedal_dll', 'onedal_core_dll']
  oneDAL_dep += cxx.find_library(lib, dirs: [pip_lib_dir])
endforeach

#add libs for tbb build dependency
tbb_dep = []
foreach lib : ['tbb12']
  tbb_dep += cxx.find_library(lib, dirs: [pip_lib_dir])
endforeach

#add libs for mkl build dependency
mkl_dep = []
foreach lib : ['mkl_intel_ilp64_dll', 'mkl_tbb_thread_dll', 'mkl_core_dll']
  mkl_dep += cxx.find_library(lib, dirs: [pip_lib_dir])
endforeach

#add link args
link_args = []
if WIN
  link_args = [
    link_args,
    '-lws2_32',   #ws2_32.lib from Windows SDK 
    '-lcrypt32',  #crypt32.lib from Windows SDK
    '-lsecur32',  #secur32.lib from Windows SDK
    '-lwldap32',  #wldap32.lib from Windows SDK
    '-lwbemuuid'  #wbemuuid.lib from Windows SDK
  ]
endif
add_global_link_arguments('-L' + lib_dir, language: 'cpp')

dlls = [
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

  #tbb dlls
  join_paths(pip_bin_dir, 'tbb12.dll'),

  #mkl dlls
  join_paths(pip_bin_dir, 'mkl_core.2.dll'),
  join_paths(pip_bin_dir, 'mkl_tbb_thread.2.dll'),
  join_paths(pip_bin_dir, 'mkl_blacs_intelmpi_ilp64.2.dll')
]

foreach dll : dlls
  install_data(dll)
endforeach

subdir('sgs')
