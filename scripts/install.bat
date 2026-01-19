::update git submodule
git submodule update --init --recursive

:: run vcpkg initialization commands
cd sgspy\extern\vcpkg
Powershell -Command .\bootstrap-vcpkg
.\vcpkg install boost-asio:x64-windows-release
.\vcpkg install gdal:x64-windows-release
cd ../../..

