#!/bin/bash
# Exit on error.
set -e

# ------------------------------------------------------------------------------------------------ #
# Default build mode is release, if "export VESTEC_DEBUG_BUILD=true" is executed before, the       #
# application will be built in debug mode.                                                         #
# Usage:                                                                                           #
#    ./make.sh [additional CMake flags, defaults to -G "Eclipse CDT4 - Unix Makefiles"]            #
# Examples:                                                                                        #
#    ./make.sh                                                                                     #
#    ./make.sh -G "Unix Makefiles"                                                                 #
#    ./make.sh -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \           #
#                      -DCMAKE_C_COMPILER_LAUNCHER=ccache                                          #
# ------------------------------------------------------------------------------------------------ #

# create some required variables -------------------------------------------------------------------

# The CMake generator and other flags can be passed as parameters.
CMAKE_FLAGS=(-G "Eclipse CDT4 - Unix Makefiles")
if [ $# -ne 0 ]; then
  CMAKE_FLAGS=( "$@" )
fi

# Check if VESTEC debug build is enabled with "export VESTEC_DEBUG_BUILD=true".
BUILD_TYPE=release
case "$VESTEC_DEBUG_BUILD" in
  (true) echo "VESTEC debug build is enabled!"; BUILD_TYPE=debug;
esac

# This directory should contain the top-level CMakeLists.txt - it is assumed to reside in the same
# directory as this script.
CMAKE_DIR="$( cd "$( dirname "$0" )" && pwd )"

# Get the current directory - this is the default location for the build and install directory.
CURRENT_DIR="$(pwd)"

# The build directory.
BUILD_DIR="$CURRENT_DIR/build/linux-$BUILD_TYPE"

# The install directory.
INSTALL_DIR="$CURRENT_DIR/install/linux-$BUILD_TYPE"

# This directory should be used as the install directory for make_externals.sh.
EXTERNALS_INSTALL_DIR="$CURRENT_DIR/install/linux-externals-$BUILD_TYPE"

# create the build directory if necessary ----------------------------------------------------------

if [ ! -d "$BUILD_DIR" ]; then
  mkdir -p "$BUILD_DIR"
fi

# configure, compile & install ---------------------------------------------------------------------
export Paraview_DIR=$EXTERNALS_INSTALL_DIR/lib/cmake/paraview-5.6

cd "$BUILD_DIR"
cmake "${CMAKE_FLAGS[@]}" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCOSMOSCOUT_EXTERNALS_DIR="$EXTERNALS_INSTALL_DIR" \
      -DParaView_DIR=$Paraview_DIR \
      -DUSE_CATALYST=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=On "$CMAKE_DIR"

cmake --build . --target install --parallel 8

# Delete empty files installed by cmake
find "$INSTALL_DIR" -type d -empty -delete

echo "Finished successfully."
