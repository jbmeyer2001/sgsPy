#update git submodule
git submodule update --init --recursive

#run vcpkg intialization commands
cd sgspy/extern/vcpkg
./bootstrap-vcpkg.sh
./vcpkg install pkgconf
./vcpkg install boost-asio
./vcpkg install gdal
cd ../../..

