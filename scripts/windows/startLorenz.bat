@echo off
set PATH=../lib:%PATH%
set PYTHONPATH=%PYTHONPATH%:Lib/site-packages
set PYTHONPATH=%PYTHONPATH%:../bin:../lib

set OMP_NUM_THREADS=%2
set OMP_DISPLAY_ENV=true
set OMP_PROC_BIND=spread #close,master,spread
set OMP_PLACES=cores #socket,cores,threads

start smpd -d 3
mpiexec -hosts 1 localhost -cores %1 VestecCatalystEmulator.exe "C:\VESTEC\Datasets\lorenz" .vtk 1 criticalPoints_Lorenz.py
@echo on
