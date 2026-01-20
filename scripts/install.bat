::update git submodule
git submodule update --init --recursive

:: run vcpkg initialization commands
cd sgspy\extern\vcpkg
Powershell -Command .\bootstrap-vcpkg
.\vcpkg install boost-asio
.\vcpkg install gdal
.\vcpkg install pkgconf
cd ../../..

