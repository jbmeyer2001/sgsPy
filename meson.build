project(
  'sgs',
  'cpp',
  default_options: [
    'cpp_std=c++20'
  ],
)

#get build operating system
WIN = host_machine.system() == 'windows' or host_machine.system() == 'cygwin'

#ensure intel oneDAL git submodule is updated
run_command('git', 'submodule', 'update', '--init', '--recursive')

#install required dependencies using vcpkg
if WIN
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.bat', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'boost-asio:x64-windows-static', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11:x64-windows-static', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'intel-mkl:x64-windows', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'tbb', check: true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal:x64-windows-static', '--recurse', check:true)
else #linux
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'boost-asio', check: true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11', check:true)
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
else #linux
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-linux/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib'
  intel_lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib/intel64'
endif
oneDAL_inc_dir = include_directories('sgs/extern/oneDAL/cpp/daal/include')
utils_inc_dir = include_directories('sgs/utils')
balanced_sampling_inc_dir = include_directories('sgs/extern/BalancedSampling')

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

  tbb_libs = [
    'tbb12',
    'tbbbind_2_5',
    'tbbmalloc',
    'tbbmalloc_proxy',
  ]

  foreach lib : mkl_libs
    lib_deps += cxx.find_library(lib, dirs : [mkl_lib_dir])
  endforeach

  foreach lib : tbb_libs
    lib_deps += cxx.find_library(lib, dirs : [tbb_lib_dir])
  endforeach

else #linux
  python_lib = run_command('find', lib_dir, '-regex', '.*libpython.*a', check:true).stdout().strip().split('/').get(-1).replace('lib', '').replace('.a', '')

  #in linux list all of the libs rather than programatically finding them
  #because of multiple definition errors. 
  libs = [
    python_lib,
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
    'mkl_blacs_intelmpi_lp64',
    'mkl_blas95_lp64',
    'mkl_cdft_core',
    'mkl_core',
    'mkl_intel_lp64',
    'mkl_lapack95_lp64',
    'mkl_scalapack_lp64',
    'mkl_tbb_thread',
    'tbb',
    'tbbmalloc'
  ]
  
  foreach lib : libs
    lib_deps += cxx.find_library(lib, dirs : [lib_dir, intel_lib_dir], static:true)
  endforeach
endif

py = import('python').find_installation(pure: false)

subdir('sgs')
