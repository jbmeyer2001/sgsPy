project('sgs', 'cpp')

run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'pybind11', check:true)

gdal_dep = dependency('gdal', required : true)
inc_dir = include_directories('sgs/extern/vcpkg/installed/x64-linux/include')

py = import('python').find_installation(pure: false)

subdir('sgs')
