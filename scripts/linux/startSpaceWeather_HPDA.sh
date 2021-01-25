#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export PV_PLUGIN_PATH=../lib/paraview-5.8/plugins/VestecPlugins

export OMP_NUM_THREADS=$3
export OMP_DISPLAY_ENV=true
export OMP_PROC_BIND=spread #close,master,spread
export OMP_PLACES=cores #socket,cores,threads

srun -N $1 -n $2 -c 56 --exclusive ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/space_weather/run007/ .vtk 0.02 $PWD/scripts/criticalPoints_SpaceWeather.py 
#> criticalPoints_spaceWeather_$1_$2_$3.csv


