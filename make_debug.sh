#!/bin/bash

# ------------------------------------------------------------------------------
# usage: ./make_release.sh [build_directory] [install_directory]
#
# This script loads the required modules, executes cmake, make and make install.
# ------------------------------------------------------------------------------

# this directory should contain the top-level CMakeLists.txt - it is assumed to
# reside in the same directory as this script
CMAKE_DIR="$( cd "$( dirname "$0" )" && pwd )"

# get the current directory - this is the default location for the build and
# install directory
CURRENT_DIR=$(pwd)

# the build directory can be passed as first parameter
BUILD_DIR="${1:-$CURRENT_DIR/build/linux-debug}"

# the install directory can be passed as second parameter
INSTALL_DIR="${2:-$CURRENT_DIR/install/linux-debug}"

# This directory should be the one used as install directory for make_externals.sh.
EXTERNALS_INSTALL_DIR="${3:-$CURRENT_DIR/install/linux-externals}"

# create build directory if neccessary -----------------------------------------

if [ ! -d $BUILD_DIR ]; then
  mkdir -p $BUILD_DIR
fi

# configure cmake --------------------------------------------------------------
export Paraview_DIR=$EXTERNALS_INSTALL_DIR/lib/cmake/paraview-5.6

cd $BUILD_DIR
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug \
-DCMAKE_EXPORT_COMPILE_COMMANDS=YES \
-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
-DEXTERNALS_DIR=$EXTERNALS_INSTALL_DIR \
-DParaView_DIR=$Paraview_DIR \
-DUSE_CATALYST=ON \
$CMAKE_DIR 


# compilation & installation ---------------------------------------------------
make -j install

