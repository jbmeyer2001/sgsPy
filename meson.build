project(
  'sgs',
  'cpp',
  default_options: [
    'cpp_std=c++17'
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
pip_build_dir = run_command(py, '-c', 'import os; print(os.environ["PATH"])', check: true).stdout().strip().split(':').get(0).split('bin').get(0)
pip_inc_dir = pip_build_dir + 'include'
#oneapi_inc_dir = pip_build_dir + 'include/oneapi'
#dal_inc_dir = pip_build_dir + 'include/dal'
pip_lib_dir = pip_build_dir + 'lib'

run_command('git', 'submodule', 'update', '--init', '--recursive')

#install required dependencies using vcpkg
if WIN
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.bat', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'boost-asio:x64-windows-static', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'intel-mkl', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'tbb', check: true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal:x64-windows-static', '--recurse', check:true)
else #linux
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'boost-asio', check: true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'intel-mkl', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'tbb', check: true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal', '--recurse', check:true)
endif

#define directories
if WIN
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-windows-static/include')
  mkl_inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-windows/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-windows-static/lib'
  mkl_lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-windows/lib/intel64'
  tbb_lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-windows/lib'
  mkl_tbb_bin = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-windows/bin' 
else #linux
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-linux/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib'
  intel_lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib/intel64'
endif
utils_inc_dir = include_directories('sgs/utils')
#oneapi_inc_dir = include_directories('sgs/extern/oneDAL/cpp')
#oneapi_daal_inc_dir = include_directories('sgs/extern/oneDAL/cpp/daal')
#oneapi_daal_include_inc_dir = include_directories('sgs/extern/oneDAL/cpp/daal/include')
root_dir = meson.current_source_dir()

#add to link args if necessary
if WIN
  win_link_args = [
    '-lws2_32',   #ws2_32.lib from Windows SDK 
    '-lcrypt32',  #crypt32.lib from Windows SDK
    '-lsecur32',  #secur32.lib from Windows SDK
    '-lwldap32',  #wldap32.lib from Windows SDK
    '-lwbemuuid'  #wbemuuid.lib from Windows SDK
  ]
endif

#get compiler, and use to find lib files in lib directory
cxx = meson.get_compiler('cpp')

lib_deps = []
mkl_lib_deps = []
oneDAL_lib_deps = []
tbb_lib_deps = []

#add all libraries in the vcpkg lib directory to the project
#the following is hacky meson build code to add all the external libraries required for the project from the vcpkg install folder
if WIN
  #the 'tree' command *should* work with both PowerShell and CMD
  libs = run_command('tree', 'sgs/extern/vcpkg/installed/x64-windows-static/lib', '/f', '/a', check: true).stdout().strip().split('|')
  foreach lib : libs
    lib = lib.strip()

    #ensure the file is actually a library, since the tree command also searches subdirectories
    #and there are often files of other types which we don't want
    if lib.endswith('.lib')

      #remove 'lib' at beginning of library if it exists
      if lib.startswith('lib')
        lib = lib.replace('lib', '')
      endif

      #get just the lib name and add it to the dependencies
      lib = lib.split('.').get(0)
      lib_deps += cxx.find_library(lib, dirs: lib_dir, static: true)
   
    endif
  endforeach

  mkl_libs = [
    'mkl_blas95_ilp64',
    'mkl_core_dll',
    'mkl_intel_ilp64_dll',
    'mkl_intel_thread_dll',
    'mkl_lapack95_ilp64',
  ]

  #tbb_libs = [
  #    ]

  #dlls = [
  #  'mkl_core.2',
  #  'mkl_tbb_threadd.2',
  #  'tbb12',
  #  'tbbbind_2_5',
  #  'tbbmalloc',
  #  'tbbmalloc_proxy',
  #]

  foreach lib : mkl_libs
    lib_deps += cxx.find_library(lib, dirs : [mkl_lib_dir], static:true)
  endforeach

  #foreach lib : tbb_libs
  #  lib_deps += cxx.find_library(lib, dirs : [tbb_lib_dir])
  #endforeach

  #foreach dll : dlls
  #  lib_deps += cxx.find_library(lib, dirs : [mkl_tbb_bin])
  #endforeach

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


  #oneDAL_libs_TEST = [
    #static_library('libonedal.a', pip_lib_dir + '/libonedal.a', implicit_include_directories: false),
    #static_library('libonedal_core.a', pip_lib_dir + '/libonedal_core.a'),
    #static_library('libonedal_parameters.a', pip_lib_dir + '/libonedal_parameters.a'),
    #static_library('libonedal_thread.a', pip_lib_dir + '/libonedal_thread.a'),
    #shared_library('libhwloc.so', pip_lib_dir + '/libhwloc.so'),
    #shared_library('libhwloc.so.15', pip_lib_dir + '/libhwloc.so.15'),
    #shared_library('libhwloc.so.15.8.1', pip_lib_dir + '/libhwloc.so.15.8.1'),
    #shared_library('libtbb.so', pip_lib_dir + '/libtbb.so'),
    #shared_library('libtbb.so.12', pip_lib_dir + '/libtbb.so.12'),
    #shared_library('libtbb.so.12.16', pip_lib_dir + '/libtbb.so.12.16'),
    #shared_library('libtbbbind_2_5.so', pip_lib_dir + '/libtbbbind_2_5.so'),
    #shared_library('libtbbbind_2_5.so.3', pip_lib_dir + '/libtbbbind_2_5.so.3'),
    #shared_library('libtbbbind_2_5.so.3.16', pip_lib_dir + '/libtbbbind_2_5.so.3.16'),
    #shared_library('libtbbmalloc.so', pip_lib_dir + '/libtbbmalloc.so'),
    #shared_library('libtbbmalloc.so.2', pip_lib_dir + '/libtbbmalloc.so.2'),
    #shared_library('libtbbmalloc.so.2.16', pip_lib_dir + '/libtbbmalloc.so.2.16'),
    #shared_library('libtbbmalloc_proxy.so', pip_lib_dir + '/libtbbmalloc_proxy.so'),
    #shared_library('libtbbmalloc_proxy.so.2', pip_lib_dir + '/libtbbmalloc_proxy.so.2'),
    #shared_library('libtbbmalloc_proxy.so.2.16', pip_lib_dir + '/libtbbmalloc_proxy.so.2.16'),
    #shared_library('liblibtcm.so', pip_lib_dir + '/libtcm.so'),
    #shared_library('libtcm.so.1', pip_lib_dir + '/libtcm.so.1'),
    #shared_library('libtcm.so.1.4.0', pip_lib_dir + '/libtcm.so.1.4.0')
    #]

  oneDAL_libs = [
    'onedal',
    'onedal_core',
    'onedal_parameters',
    'onedal_thread'
  ]

  #message('py install dir:')
  #message(py.get_install_dir())
  
  foreach lib : oneDAL_libs
    oneDAL_lib_deps += cxx.find_library(lib, dirs: [pip_lib_dir])
  endforeach
  
  #install_data(pip_lib_dir + '/libtbb.so.12', install_dir: py.get_install_dir())
  
  #oneDAL_dep = declare_dependency(
  #  include_directories: pip_inc_dir,
  #  link_with: oneDAL_libs
  #)

endif

pybind11_dep = dependency('pybind11')

subdir('sgs')
