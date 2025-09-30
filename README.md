## Overview
The sgsPy Python pacakge (Structurally Guided Sampling) is a package in development.

While the package has not yet been officially distributed or released, beta testers
are highly encouraged! See the lower sections of this readme for installation and
running instructions. 

#### How to stay up to date?
Either by watching this repository, or if you would like to recieve more detailed information on
larger updates, you may email me at jmeyer03 'at' mail.ubc.ca asking to join the sgs mailing list.

#### How to contribute?
The easiest way to contribute is by logging any errors, or posting any ideas as an issue on this github repository.

#### known issues:
 - the implementation of sample.strat does not work on windows. A new method which is effective across large areas, and runs error-free on windows is under development.
 - occasionally, an access-related test will fail.

## Installation
### NOTE: when distributed, the package will be available through PyPI and thus a quick 'pip install sgs'. For now, the package can be installed and built as follows depending on operating system.

### Linux:
1. Ensure you have Python with pip installed, and git.
2. clone the repository with the following command:
```
git clone https://github.com/jbmeyer2001/sgsPy.git
```

3. create a python virtual environment (highly recommended)
```
python -m venv .venv
```

the virtual environment then resides in the .venv folder, and can be activated with the following command:
```
source ./.venv/bin/activate
```

4. install dependency requirements. If you do not already have these, the dependency installation will fail.
 - build-essential (if you don't already have a C++ compiler)
 - pkg-config
 - auto-conf
 - libtool
 - bison
 - flex

```
sudo apt install build-essential
sudo apt install pkg-config
sudo apt install auto-conf
sudo apt install libtool
sudo apt install bison
sudo apt install flex
```

5. Run the following commands to install dependencies from vcpkg (a C++ package manager). These commands will be automatically run any time the project is built, but the first time it is recommended to run them seperately to view the outputs. These commands may take a while (potentially a few hours), but once they are installed the build command will run much faster.
```
git submodule update --init --recursive #check out vcpkg submodule
cd sgs/extern/vcpkg
./bootstrap-vcpkg
./vcpkg integrate install
./vcpkg install gdal
./vcpkg install boost-asio
./vpckg install intel-mkl
./vcpkg install tbb
```

6. build the package: run the following command from within the folder containing this file.
```
pip install .
```

7. If you intend to run the tests (using pytest), both pytest and geopandas are required and can be installed as follow:
```
pip install pytest
pip install geopandas
```

### Windows:

1. ensure you have Python with pip, git, and a C++ compiler.

2. If you do not already have a C++ compiler, one can installed by installing the Microsoft Visual Studio IDE with C++ build tools.

3. clone the repository with the following command.
```
git clone https://github.com/jbmeyer2001/sgsPy.git
```

4. create a Python virtual environment (highly recommended).
```
python -m venv .venv
```

if you are usign powershell, the activation command is shown below, the activation command for command prompt may be different.
```
./.venv/Scripts/activate
```

5. Run the following commands to install dependencies from vcpkg (a C++ package manager). These commands will be automatically run any time the project is built, but the first time it is recommended to run them seperately to view the outputs. These commands may take a while (potentially a few hours), but once they are installed the build command will run much faster.
```
git submodule update --init --recursive #check out vcpkg submodule
cd sgs/extern/vcpkg
./bootstrap-vcpkg
./vcpkg integrate install
./vcpkg install gdal
./vcpkg install boost-asio
./vpckg install intel-mkl
./vcpkg install tbb
```

6. build the package: run the following command from within the folder containing this file.
```
pip install .
```

7. If you intend to run the tests (using pytest), both pytest and geopandas are required and can be installed as follow:
```
pip install pytest
pip install geopandas
```

## How to run sgsPy

Tests may be ran by running the following command from within the folder containing this file. Both pytest and geopandas must be installed to run the tests.
```
pytest
```

To run sgs on your own data, the best resource to look to is the tests, which are in the 'tests/' directory. Additionally, preliminary documentation for each function is available within the Python file that contains that function. For example, for information on how to run quantiles stratification on a raster, look to the comments within 'sgs/stratify/quantiles/quantiles.py'. Alternatively, you may contact me directly at jmeyer03 'at' mail.ubc.ca . 

## Develoment progress
### stratification:
 - breaks (including large raster processing)
 - quantiles (including large raster processing)
 - map (including large raster processing)
 - poly (including large raster processing)

### sampling:
 - stratified sampling (NOT including large raster processing)
 - simple random sampling (NOT including large raster processing)
 - systematic sampling (NOT including large raster processing)

### Projected future development priorities:
 - upgrade sampling methods to process large rasters
 - implement principal component analysis for large rasters
 - add existing sample plots as parameter for sampling methods
 - debug occasionally failing access-related tests
 - implement kmeans stratification

## Licensing
the sgsPy package is licensed under the MIT license. the LICENSE.txt file contains
full licensing information for all dependencies. None of the licenses are
copyleft licenses sush as GPL. However, it's worth noting that two dependencies,
oneMKL and oneTBB, are licensed under the Intel Simplified Software License (ISSL).
The ISSL allows use and redistribution, but is not technically open source since the
code is propriatery. As such, this software cannot be distributed alongside other
software which uses a copyleft license such as GPL, since oneMKL and oneTBB are not
technically open source.
