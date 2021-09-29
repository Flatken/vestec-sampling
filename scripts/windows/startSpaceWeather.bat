@echo off

set DATA_PATH="D:\\vr_data\\VESTEC\\space_weather\\run007"
set PV_PLUGIN_PATH=.;paraview-5.9/plugins/VestecPlugins;paraview-5.9/plugins/TopologyToolKit
set NUM_THREADS=6
set MPI_PROCESSES=1

start smpd -d 3
mpiexec -env OMP_NUM_THREADS %NUM_THREADS% -env PV_PLUGIN_PATH %PV_PLUGIN_PATH% -hosts 1 localhost -cores %MPI_PROCESSES% VestecCatalystEmulator.exe %DATA_PATH% .vtk 12 scripts\criticalPoints_SpaceWeather.py
@echo on
