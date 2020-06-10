@echo off
set PATH=../lib:%PATH%
set PYTHONPATH=%PYTHONPATH%:Lib/site-packages
set PYTHONPATH=%PYTHONPATH%:../bin:../lib
start smpd -d 3
mpiexec -hosts 1 localhost -cores 1 VestecCatalystEmulator.exe "D:\\vr_data\\TaylorGreen\\" .vtk 1 criticalPoints_Evaluation.py
@echo on
