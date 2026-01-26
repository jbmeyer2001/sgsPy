#update git submodule
git submodule update --init --recursive

#set pkgconfig environment variables
export PKG_CONFIG="${PWD}/sgspy/extern/vcpkg/installed/x64-linux/tools/pkgconf/pkgconf"
export PKG_CONFIG_PATH="${PWD}/sgspy/extern/vcpkg/installed/x64-linux/lib/pkgconfig"

#run build command -- requires the 'build' python package downloadable with 'pip install build'
clear
python -m build --wheel --outdir dist 

#run auditwheel on wheel to ensure manylinux compliance, excluding libraries which will come from a Python dependency
auditwheel repair dist/*.whl -w dist --plat manylinux_2_35_x86_64 --exclude libmkl_core.so.2 --exclude libmkl_tbb_thread.so.2 --exclude libmkl_intel_ilp64.so.2 --exclude libonedal_parameters.so.3 --exclude libonedal_core.so.3 --exclude libonedal_thread.so.3 --exclude libonedal.so.3
rm dist/*linux_x86_64.whl
