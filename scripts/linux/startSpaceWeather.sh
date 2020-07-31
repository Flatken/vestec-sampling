#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
mpirun.mpich -np 1 ./VestecCatalystEmulator /unsecured/flat_ma/vestec/datasets/space_weather/run007/ .vtk 0.02 criticalPoints_SpaceWeather.py

