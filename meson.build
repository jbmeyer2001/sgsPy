project(
  'sgs', 
  'cpp',
  default_options: [
      'b_vscrt=mt'
  ],
)

#get build operating system
if host_machine.system() == 'windows' or host_machine.system() == 'cygwin'
  WIN = true
endif

#install required dependencies using vcpkg
if WIN
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.bat', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11:x64-windows-static', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal:x64-windows-static', '--recurse', check:true)
else
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal', '--recurse', check:true)
endif

#define directories
if WIN
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-windows-static/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-windows-static/lib'
else
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-linux/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib'
endif
utils_inc_dir = include_directories('sgs/utils')
balanced_sampling_inc_dir = include_directories('sgs/extern/balancedSampling')

#get the python lib name for any version of python
if WIN
  #the following is hacky meson build code to get a version agnostic python library name, run from either powershell or cmd
  tree_output = run_command('tree', 'sgs/extern/vcpkg/installed/x64-windows-static/lib', '/f', '/a', check: true).stdout().strip().split('|')
  foreach filename : tree_output
    stripped_filename = filename.strip()
    if stripped_filename.startswith('python') and stripped_filename.endswith('.lib')
      python_lib = stripped_filename.split('.').get(0)
    endif
  endforeach
else
  python_lib = run_command('find', lib_dir, '-regex', '.*libpython.*a', check:true).stdout().strip().split('/').get(-1).replace('lib', '').replace('.a', '')
endif

win_libs = [
  'aec-static',
  'charset',
  'iconv',
  'expatMT',
  'qhullstatic',
  'szip-static',
  'zlib', 
]

linux_libs = [
  'aec',
  'expat',
  'pkgconf',
  'png',
  'sz',
  'z',
  'mod_spatialite',
]

common_libs = [
  python_lib,
  'crypto',
  'curl',
  'ecpg',
  'ecpg_compat',
  'hdf5',
  'hdf5_cpp',
  'hdf5_hl',
  'hdf5_hl_cpp',
  'pgcommon',
  'pgport',
  'pgtypes',
  'png16',
  'pq',
  'ssl',
  'sharpyuv',
  'webp',
  'webpdecoder',
  'webpdemux',
  'webpmux',
  'xml2',
  'freexl',
  'gdal',
  'geos',
  'geos_c',
  'geotiff',
  'gif',
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
  'proj',
  'qhullcpp',
  'qhullstatic_r',
  'spatialite',
  'sqlite3',
  'tiff',
  'turbojpeg',
  'uriparser',
  'zstd',
]

lib_deps = []

#get compiler, and use to find lib files in lib directory
cxx = meson.get_compiler('cpp')
foreach lib : common_libs
  lib_deps += cxx.find_library(lib, dirs : lib_dir, static:true)
endforeach

if WIN
  foreach lib : win_libs
    lib_deps += cxx.find_library(lib, dirs : lib_dir, static:true)
  endforeach
else
  foreach lib : linux_libs
    lib_deps += cxx.find_library(lib, dirs : lib_dir, static:true)
  endforeach
endif

py = import('python').find_installation(pure: false)

subdir('sgs')
