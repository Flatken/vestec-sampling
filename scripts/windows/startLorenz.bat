@echo off

set DATA_PATH="C:\\VESTEC\\Datasets\\lorenz"
@REM set DATA_PATH="D:\\vr_data\\LorenzAttractor"

set PV_PLUGIN_PATH=.;paraview-5.9/plugins/VestecPlugins;paraview-5.9/plugins/TopologyToolKit
set NUM_THREADS=6
set MPI_PROCESSES=1

start smpd -d 3
mpiexec -env PV_PLUGIN_PATH %PV_PLUGIN_PATH% -env OMP_NUM_THREADS %NUM_THREADS% -env OMP_DISPLAY_ENV true -env OMP_PROC_BIND spread -env OMP_PLACES cores -hosts 1 localhost -cores %MPI_PROCESSES% VestecCatalystEmulator.exe %DATA_PATH% .vtk 1 scripts\criticalPoints_Lorenz.py
@echo on
