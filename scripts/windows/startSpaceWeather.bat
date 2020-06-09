@echo off
set PATH=../lib:%PATH%
set PYTHONPATH=%PYTHONPATH%:Lib/site-packages
set PYTHONPATH=%PYTHONPATH%:../bin:../lib
start smpd -d 3
mpiexec -hosts 1 localhost -cores 1 VestecCatalystEmulator.exe "C:\\VESTEC\\Datasets\\space_weather\\run007" .vtk 1 criticalPoints_SpaceWeather.py
@echo on
