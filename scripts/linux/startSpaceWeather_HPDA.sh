#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$3
srun -N $1 -n $2 -c 56 --exclusive ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/space_weather/run007/ .vtk 0.02 criticalPoints_SpaceWeather.py > criticalPoints_spaceWeather_$1_$2_$3.csv


