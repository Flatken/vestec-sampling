#!/bin/bash

# exit on error
set -e

# Create some required variables. ------------------------------------------------------------------

# This directory should all the submodules - they are assumed to reside in the subdirectory 
# "externals" next to this script.
EXTERNALS_DIR="$( cd "$( dirname "$0" )" && pwd )/externals"

# Get the current directory - this is the default location for the build and install directory.
CURRENT_DIR="$(pwd)"

# Check if VESTEC debug build is enabled with "export VESTEC_DEBUG_BUILD=true".
BUILD_TYPE=release
case "$VESTEC_DEBUG_BUILD" in
  (true) echo "VESTEC debug build is enabled!"; BUILD_TYPE=debug;
esac

# The build directory can be passed as first parameter.
BUILD_DIR="${1:-$CURRENT_DIR/build/linux-externals-$BUILD_TYPE}"

# The install directory can be passed as second parameter.
INSTALL_DIR="${2:-$CURRENT_DIR/install/linux-externals-$BUILD_TYPE}"

# Create some default installation directories.
cmake -E make_directory "$INSTALL_DIR/lib"
cmake -E make_directory "$INSTALL_DIR/share"
cmake -E make_directory "$INSTALL_DIR/bin"
cmake -E make_directory "$INSTALL_DIR/include"

# MPIR -----------------------------------------------------------------------------------------
echo ""
echo "Building and installing MPIR Library ..."

mkdir $BUILD_DIR/mpir
cd $EXTERNALS_DIR/mpir
bash autoconf.sh
cd $BUILD_DIR/mpir
$EXTERNALS_DIR/mpir/configure -prefix $INSTALL_DIR
make
make install

echo ""

cd "$CURRENT_DIR"

## --------------------------------------------------------------------------------------------------
echo "Finished successfully."
#
