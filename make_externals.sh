#!/bin/bash

# exit on error
set -e

# Create some required variables. ------------------------------------------------------------------

# This directory should all the submodules - they are assumed to reside in the subdirectory 
# "externals" next to this script.
EXTERNALS_DIR="$( cd "$( dirname "$0" )" && pwd )/externals"

# Get the current directory - this is the default location for the build and install directory.
CURRENT_DIR="$(pwd)"

# The build directory can be passed as first parameter.
BUILD_DIR="${1:-$CURRENT_DIR/build/linux-externals}"

# The install directory can be passed as second parameter.
INSTALL_DIR="${2:-$CURRENT_DIR/install/linux-externals}"

# Create some default installation directories.
cmake -E make_directory "$INSTALL_DIR/lib"
cmake -E make_directory "$INSTALL_DIR/share"
cmake -E make_directory "$INSTALL_DIR/bin"
cmake -E make_directory "$INSTALL_DIR/include"

# TTK patch paraview--------------------------------------------------------------------------------
#cd $EXTERNALS_DIR/ttk/paraview/patch/
#./patch-paraview-5.6.0.sh $EXTERNALS_DIR/paraview-5.6

# Paraview -----------------------------------------------------------------------------------------
echo "Building and installing Paraview 5.6 ..."
echo ""
echo ""

cmake -E make_directory "$BUILD_DIR/paraview" && cd "$BUILD_DIR/paraview"
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DPARAVIEW_ENABLE_PYTHON=ON \
      -DVTK_PYTHON_VERSION=3 \
      -DPARAVIEW_BUILD_QT_GUI=On \
      -DPARAVIEW_USE_MPI=ON \
      -DPARAVIEW_ENABLE_CATALYST=ON \
      -DCMAKE_BUILD_TYPE=Release "$EXTERNALS_DIR/paraview-5.6" 
cmake --build . --target install --parallel 8

# TTK -------------------------------------------------------------------------------------------

echo ""
echo "Building and installing TTK ..."
echo ""

cmake -E make_directory "$BUILD_DIR/ttk" && cd "$BUILD_DIR/ttk"
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_PREFIX_PATH=$INSTALL_DIR/lib/cmake/paraview-5.6 \
      -DParaView_CMAKE_DIR=$INSTALL_DIR/lib/cmake/paraview-5.6 \
      -DCMAKE_BUILD_TYPE=Release "$EXTERNALS_DIR/ttk"
cmake --build . --target install --parallel 8

## --------------------------------------------------------------------------------------------------
echo "Finished successfully."
#