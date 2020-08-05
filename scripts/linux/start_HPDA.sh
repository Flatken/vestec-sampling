#!/bin/bash

export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH
srun -N $1 -c $2 --exclusive ./VestecCatalystEmulator /scratch/VESTEC-DATASETS/2D/ .vtk 0.02 criticalPoints_Evaluation.py > $3

