#!/bin/sh

# This is a simple script that runs the minimum commands needed to 
# compile and install Sire, assuming it has already been successfully
# compiled and installed. This is useful to speed up development, 
# as you can run this script instead of `compile_sire.sh` in the 
# above directory.

# Only run this script if you know what you are doing ;-)

set -e

cd corelib
$HOME/sire.app/bin/cmake --build . --target install -- VERBOSE=1 -j 8
cd ../wrapper
$HOME/sire.app/bin/cmake --build . --target install -- VERBOSE=1 -j 8
cd ..
