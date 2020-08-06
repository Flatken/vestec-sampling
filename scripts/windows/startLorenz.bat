@echo off
set PATH=../lib:%PATH%
set PYTHONPATH=%PYTHONPATH%:Lib/site-packages
set PYTHONPATH=%PYTHONPATH%:../bin:../lib
start smpd -d 3
mpiexec -hosts 1 localhost -cores 4 VestecCatalystEmulator.exe "D:\\vr_data\\LorenzAttractor\\" .vtk 1 criticalPoints_Lorenz.py
@echo on
