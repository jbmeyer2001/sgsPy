Statement of ongoing development
The sgsPy Python pacakge (Structurally Guided Sampling) is a package in development.

It has not been distributed, nor officially released besides making this github 
repository public. As such, errors are not unlikely, despite any tests which it 
may pass. However, you are highly encouraged to test sgsPy out, and contribute 
to development of the package by notifying me of any bugs which you may find or
ideas you may have. The easiest way to reach me is by posting within the github
issues or by emailing jmeyer03 'at mail.ubc.ca .

Installation
clone the repository
'''https://github.com/jbmeyer2001/sgs.git'''

create a Python virtual environment (recommended)

to install dependencies and build, run the following command from the directory containing this file
'''pip install .'''

Develoment progress
stratification:
 - breaks (including large raster processing)
 - quantiles (including large raster processing)
 - map (including large raster processing)
 - poly (including large raster processing)

sampling:
 - stratified sampling (NOT including large raster processing)
 - simple random sampling (NOT including large raster processing)
 - systematic sampling (NOT including large raster processing)

Projected future development priorities:
 - upgrade sampling methods to process large rasters
 - implement principal component analysis for large rasters
 - add existing sample plots as parameter for sampling methods
 - debug occasionally failing access-related tests
 - implement kmeans stratification

Licensing
the sgsPy package is licensed under the MIT license. the LICENSE.txt file contains
full licensing information for all dependencies. None of the licenses are
copyleft licenses sush as GPL. However, it's worth noting that two dependencies,
oneMKL and oneTBB, are licensed under the Intel Simplified Software License (ISSL).
The ISSL allows use and redistribution, but is not technically open source since the
code is propriatery. As such, this software cannot be distributed alongside other
software which uses a copyleft license such as GPL, since oneMKL and oneTBB are not
technically open source.
