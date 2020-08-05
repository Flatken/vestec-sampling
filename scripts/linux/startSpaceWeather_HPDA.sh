#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$3
srun -N $1 -n $2 -c $3 --exclusive ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/space_weather/run007/ .vtk 0.02 criticalPoints_SpaceWeather.py


