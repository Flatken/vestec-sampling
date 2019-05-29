#!/bin/bash

module() { eval `/usr/bin/modulecmd bash $*`; }

export MODULEPATH=/tools/modulesystem/modulefiles

module purge
module load iv-group
module unload gdal

module load python3/python-3.7.2/sled12.x86_64.gcc.release
export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:../lib/python3.4/site-packages
export PYTHONPATH=$PYTHONPATH:../bin:../lib

echo $PYTHONPATH

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $SCRIPT_DIR

mpiexec -np 2 ./VestecCatalystEmulator /unsecured/flat_ma/vestec/datasets/test/ 0.02 $1

