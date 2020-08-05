#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
srun -N $1 -c $2 --exclusive ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/space_weather/run007/ .vtk 0.02 criticalPoints_SpaceWeather.py > $3

