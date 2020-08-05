#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
mpirun.mpich -np 1 ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/2D/ .vtk 0.02 criticalPoints_Evaluation.py

