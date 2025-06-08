project('sgs', 'cpp')

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

#add to link args if necessary
if WIN
  link_arguments = [
    '-lws2_32', 
    '-lcrypt32', 
    '-lsecur32', 
    '-lwldap32', 
    '-lwbemuuid'
  ]
else
  link_arguments = []
endif

#get compiler, and use to find lib files in lib directory
cxx = meson.get_compiler('cpp')

lib_deps = []

#add all libraries in the vcpkg lib directory to the project
if WIN
  #the following is hacky meson build code to get a version agnostic python library name, run from either powershell or cmd
  tree_output = run_command('tree', 'sgs/extern/vcpkg/installed/x64-windows-static/lib', '/f', '/a', check: true).stdout().strip().split('|')
  foreach filename : tree_output
    filename = filename.strip()
    if filename.endswith('.lib')
      if filename.startswith('lib')
        filename = filename.replace('lib', '')
      endif
      filename = filename.split('.').get(0)
      lib_deps += cxx.find_library(filename, dirs: lib_dir, static: true)
    endif
  endforeach
else
  python_lib = run_command('find', lib_dir, '-regex', '.*libpython.*a', check:true).stdout().strip().split('/').get(-1).replace('lib', '').replace('.a', '')
endif

py = import('python').find_installation(pure: false)

subdir('sgs')
