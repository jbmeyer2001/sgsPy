::update git submodule
git submodule update --init --recursive

:: edit the gdal pkgconfig file to say 'gdal.lib' instead of 'gdal'
:: this is very hacky but the build system searches for a nonexistant gdal.obj file if not done
Powershell -Command "(gc %CD%\sgspy\extern\vcpkg\installed\x64-windows\lib\pkgconfig\gdal.pc) -replace 'lib gdal', 'lib  gdal.lib' | Out-File -encoding ASCII %CD%\sgspy\extern\vcpkg\installed\x64-windows\lib\pkgconfig\gdal.pc"

:: set environment variables
set "PKG_CONFIG=%CD%\sgspy\extern\vcpkg\installed\x64-windows\tools\pkgconf\pkgconf.exe"
set "PKG_CONFIG_PATH=%CD%\sgspy\extern\vcpkg\installed\x64-windows\lib\pkgconfig"

:: run build command
cls
python -m build --wheel --outdir dist
