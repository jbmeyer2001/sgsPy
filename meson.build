project('sgs', 'cpp')

#get build operating system
WIN = host_machine.system() == 'windows' or host_machine.system() == 'cygwin'

#install required dependencies using vcpkg
if WIN
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.bat', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11:x64-windows-static', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal:x64-windows-static', '--recurse', check:true)
else #linux
  run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11', check:true)
  run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal', '--recurse', check:true)
endif

#define directories
if WIN
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-windows-static/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-windows-static/lib'
else #linux
  inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-linux/include')
  lib_dir = meson.current_source_dir() + '/sgs/extern/vcpkg/installed/x64-linux/lib'
endif
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
    'mod_spatialite'
  ]
  
  foreach lib : libs
    lib_deps += cxx.find_library(lib, dirs : [lib_dir], static:true)
  endforeach
endif

py = import('python').find_installation(pure: false)

subdir('sgs')
