#update git submodule
git submodule update --init --recursive

#run vcpkg intialization commands
cd sgspy/extern/vcpkg
./bootstrap-vcpkg.sh
./vcpkg install boost-asio:x64-linux-release
./vcpkg install gdal:x64-linux-release
cd ../../..

