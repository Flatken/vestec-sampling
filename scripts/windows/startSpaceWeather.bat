@echo off
set PATH=../lib;%PATH%
set PYTHONPATH=%PYTHONPATH%;Lib/site-packages
set PYTHONPATH=%PYTHONPATH%;../bin;../lib

start smpd -d 3
mpiexec -env OMP_NUM_THREADS %2 -env OMP_DISPLAY_ENV true -env OMP_PROC_BIND spread -env OMP_PLACES cores -hosts 1 localhost -cores %1 VestecCatalystEmulator.exe "D:\\vr_data\\VESTEC\\space_weather\\run007" .vtk 40 criticalPoints_SpaceWeather.py
@echo on
