#!/bin/bash

set -e -o pipefail

export CPATH=${PREFIX}/include
#export CMAKE_LDFLAGS="-L${PREFIX}/lib"
export CMAKE_LDFLAGS="-L${PREFIX}/lib -lXfixes -lXrender -lXext -lXft -lfontconfig -lfontconfig -lXinerama -ldl -lm -lX11"
export CPPFLAGS="$CPPFLAGS -I${PREFIX}/include"
export CXXFLAGS="-L${PREFIX}/lib -lXfixes -lXrender -lXext -lXft -lfontconfig -lfontconfig -lXinerama -ldl -lm -lX11"
export LDFLAGS="$LDFLAGS -L${PREFIX}/lib"
export LD_LIBRARY_PATH=${PREFIX}/lib
export LIBRARY_PATH=${PREFIX}/lib

mkdir -p build
cd build
# cmake -DCMAKE_INSTALL_PREFIX="${PREFIX}" -Wl -DFORCE_OWN_FFTW=ON ..
cmake -DCMAKE_INSTALL_PREFIX="${PREFIX}" -DCUDA=OFF ..
make install

