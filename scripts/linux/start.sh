#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
mpirun.mpich -np 1 ./VestecCatalystEmulator /unsecured/flat_ma/vestec/datasets/test/ .vtk 0.02 criticalPoints_Evaluation.py

