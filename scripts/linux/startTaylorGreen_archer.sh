#!/bin/bash

#SBATCH --account=d170
#SBATCH --partition=standard
#SBATCH --qos=standard 

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export PV_PLUGIN_PATH=../lib/paraview-5.8/plugins/VestecPlugins

export OMP_NUM_THREADS=$3
export OMP_DISPLAY_ENV=true
export OMP_PROC_BIND=spread #close,master,spread
export OMP_PLACES=socket #socket,cores,threads

srun -p standard -A d170 --qos=standard --time=0:40:0 -N $1 -n $2 -c 128 --exclusive ./VestecCatalystEmulator /work/d170/d170/shared/VESTEC-DATASETS/taylor_green/ .vtk 0.02 $PWD/scripts/criticalPoints_TaylorGreen.py > criticalPoints_TaylorGreen_$1_$2_$3.csv

