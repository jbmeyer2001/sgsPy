project('sgs', 'cpp')

run_command('./sgs/extern/vcpkg/bootstrap-vcpkg.sh', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'integrate', 'install', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'gdal', check:true)
run_command('./sgs/extern/vcpkg/vcpkg', 'install', 'boost', check:true)

boost_dep = dependency('boost', required : true)
gdal_dep = dependency('gdal', required : true)

py = import('python').find_installation(pure: false)

subdir('sgs')
