project('sgs', 'cpp')

run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal', '--recurse', check:true)

inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-linux/include')

cxx = meson.get_compiler('cpp')
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
  'mod_spatialite'
]
lib_deps = []

foreach lib : libs
  lib_deps += cxx.find_library(lib, dirs : [meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib'], static:true)
endforeach

py = import('python').find_installation(pure: false)

subdir('sgs')
