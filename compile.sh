#!/bin/bash
# Usage: mkdir build && cd build && cmake [options] .. && make -j 8
#
# cmake options:
#   -DCMAKE_BUILD_TYPE=Debug
#   -DUSE_NETCDF=OFF
#   -DUSE_ETA=OFF
#   -DSERIAL=ON
#   -DUSE_MPI=ON
#   -DARCH=generic|x86-64|skylake|native

export CPATH="/usr/include:/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/"
export LIBRARY_PATH="/usr/lib"
