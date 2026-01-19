:: edit the gdal pkgconfig file to say 'gdal.lib' instead of 'gdal'
:: this is very hacky but the build system searches for a nonexistant gdal.obj file if not done
Powershell -Command "(gc %CD%\sgspy\extern\vcpkg\installed\x64-windows-release\lib\pkgconfig\gdal.pc) -replace 'lib gdal', 'lib  gdal.lib' | Out-File -encoding ASCII %CD%\sgspy\extern\vcpkg\installed\x64-windows-release\lib\pkgconfig\gdal.pc"

:: set environment variables
set "PKG_CONFIG=%CD%\sgspy\extern\vcpkg\installed\x64-windows-release\tools\pkgconf\pkgconf.exe"
set "PKG_CONFIG_PATH=%CD%\sgspy\extern\vcpkg\installed\x64-windows-release\lib\pkgconfig"

:: run build command
cls
pip install .
