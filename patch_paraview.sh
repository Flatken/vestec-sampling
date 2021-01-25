#!/bin/bash

# exit on error
set -e

# Create some required variables. ------------------------------------------------------------------

# This directory should all the submodules - they are assumed to reside in the subdirectory 
# "externals" next to this script.
EXTERNALS_DIR="$( cd "$( dirname "$0" )" && pwd )/externals"

# TTK patch paraview--------------------------------------------------------------------------------
cd $EXTERNALS_DIR/ttk/paraview/patch/
./patch-paraview-5.6.0.sh $EXTERNALS_DIR/paraview-5.6

## --------------------------------------------------------------------------------------------------
echo "Finished successfully."
#
