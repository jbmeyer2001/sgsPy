#update git submodule
git submodule update --init --recursive

#set pkgconfig environment variables
export PKG_CONFIG="${PWD}/sgspy/extern/vcpkg/installed/x64-linux/tools/pkgconf/pkgconf"
export PKG_CONFIG_PATH="${PWD}/sgspy/extern/vcpkg/installed/x64-linux/lib/pkgconfig"

#run build command -- requires the 'build' python package downloadable with 'pip install build'
clear
python -m build --wheel --outdir dist 

