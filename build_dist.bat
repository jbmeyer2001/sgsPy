::update git submodule
git submodule update --init --recursive

:: run vcpkg initialization commands
cd sgs\extern\vcpkg
Powershell -Command .\bootstrap-vcpkg
.\vcpkg install boost-asio
.\vcpkg install gdal
cd ../../..

:: edit the gdal pkgconfig file to say 'gdal.lib' instead of 'gdal'
:: this is very hacky but the build system searches for a nonexistant gdal.obj file if not done
Powershell -Command "(gc %CD%\sgs\extern\vcpkg\installed\x64-windows\lib\pkgconfig\gdal.pc) -replace 'lib gdal', 'lib  gdal.lib' | Out-File -encoding ASCII %CD%\sgs\extern\vcpkg\installed\x64-windows\lib\pkgconfig\gdal.pc"

:: set environment variables
set "PKG_CONFIG=%CD%\sgs\extern\vcpkg\installed\x64-windows\tools\pkgconf\pkgconf.exe"
set "PKG_CONFIG_PATH=%CD%\sgs\extern\vcpkg\installed\x64-windows\lib\pkgconfig"

:: run build command
cls
python -m build --wheel --outdir dist
