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

# TTK patch paraview--------------------------------------------------------------------------------
#cd $EXTERNALS_DIR/ttk/paraview/patch/
#./patch-paraview-5.6.0.sh $EXTERNALS_DIR/paraview-5.6

# Paraview -----------------------------------------------------------------------------------------
echo "Building and installing Paraview ..."
echo ""
echo ""

cmake -E make_directory "$BUILD_DIR/paraview" && cd "$BUILD_DIR/paraview"
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DPARAVIEW_USE_PYTHON=ON \
      -DVTK_PYTHON_VERSION=3 \
      -DPARAVIEW_USE_QT=ON \
      -DPARAVIEW_USE_MPI=ON \
      -DPARAVIEW_USE_VTKM=OFF \
      -DCMAKE_BUILD_TYPE=Release "$EXTERNALS_DIR/paraview-5.6" 
cmake --build . --target install --parallel "$(nproc)"


# EIGEN -----------------------------------------------------------------------------------------
echo ""
echo "Building and installing Eigen ..."

cmake -E make_directory "$BUILD_DIR/eigen" && cd "$BUILD_DIR/eigen"
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DCMAKE_BUILD_TYPE=Release "$EXTERNALS_DIR/eigen"
cmake --build . --target install --parallel "$(nproc)"

echo ""

# TTK -------------------------------------------------------------------------------------------

echo ""
echo "Building and installing TTK ..."
echo ""

cmake -E remove_directory "$EXTERNALS_DIR/ttk/paraview/WRLExporter"
cmake -E remove_directory "$EXTERNALS_DIR/ttk/core/vtk/ttkWRLExporter" 

cmake -E make_directory "$BUILD_DIR/ttk" && cd "$BUILD_DIR/ttk"
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DParaView_DIR="$BUILD_DIR/paraview" \
      -DEigen3_DIR="$INSTALL_DIR/share/eigen3/cmake" \
      -DTTK_INSTALL_PLUGIN_DIR="$INSTALL_DIR/lib/plugins" \
      -DTTK_ENABLE_ZLIB=OFF \
      -DTTK_ENABLE_KAMIKAZE=On \
      -DTTK_ENABLE_MPI=ON \
      -DTTK_BUILD_STANDALONE_APPS=OFF \
      -DVTK_MODULE_ENABLE_ttkCinemaImaging=DONT_WANT \
      -DVTK_MODULE_ENABLE_ttkUserInterfaceBase=DONT_WANT \
      -DCMAKE_BUILD_TYPE=Release "$EXTERNALS_DIR/ttk"
cmake --build . --target install --parallel "$(nproc)"

echo ""

cd "$CURRENT_DIR"

## --------------------------------------------------------------------------------------------------
echo "Finished successfully."
#
