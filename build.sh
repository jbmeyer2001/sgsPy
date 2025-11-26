#update git submodule
git submodule update --init --recursive

#run vcpkg intialization commands
cd sgs/extern/vcpkg
./bootstrap-vcpkg.sh
./vcpkg install boost-asio
./vcpkg install gdal
cd ../../..

#set pkgconfig environment variables
export PKG_CONFIG="${PWD}/sgs/extern/vcpkg/installed/x64-linux/tools/pkgconf/pkgconf"
export PKG_CONFIG_PATH="${PWD}/sgs/extern/vcpkg/installed/x64-linux/lib/pkgconfig"

#run build command
clear
pip install .
