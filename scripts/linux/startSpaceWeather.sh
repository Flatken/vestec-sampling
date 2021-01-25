#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
export PV_PLUGIN_PATH=../lib/paraview-5.8/plugins/VestecPlugins
mpirun.mpich -np 1 ./VestecCatalystEmulator /unsecured/flat_ma/vestec/datasets/space_weather/run007/ .vtk 0.02 $PWD/scripts/criticalPoints_SpaceWeather.py

