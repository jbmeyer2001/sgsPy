#set pkgconfig environment variables
export PKG_CONFIG="${PWD}/sgspy/extern/vcpkg/installed/x64-linux/tools/pkgconf/pkgconf"
export PKG_CONFIG_PATH="${PWD}/sgspy/extern/vcpkg/installed/x64-linux/lib/pkgconfig"

#run build command
clear
pip install .
