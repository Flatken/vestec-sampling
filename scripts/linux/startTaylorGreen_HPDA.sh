#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export PV_PLUGIN_PATH=../lib/paraview-5.8/plugins/VestecPlugins

export OMP_NUM_THREADS=$3
export OMP_DISPLAY_ENV=true
export OMP_PROC_BIND=true #close,master,spread,true,false
export OMP_PLACES=cores #sockets,cores,threads
export OMP_SCHEDULE=static #dynamic,static

export SLURM_CPU_BIND=verbose

srun -N $1 -n $2 -c 56 --cpu-bind=no ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/taylor_green/ .vtk 0.02 $PWD/scripts/criticalPoints_TaylorGreen.py > criticalPoints_TaylorGreen_$1_$2_$3.csv


