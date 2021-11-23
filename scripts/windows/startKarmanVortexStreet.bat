@echo off

set DATA_PATH="C:\\VESTEC\\Datasets\\time_varying_flow"
@REM set DATA_PATH="D:\\vr_data\\KarmanVortexStreet\\2D\\"
set PV_PLUGIN_PATH=.;paraview-5.9/plugins/VestecPlugins;paraview-5.9/plugins/TopologyToolKit
set NUM_THREADS=6
set MPI_PROCESSES=1

start smpd -d 3
mpiexec -env OMP_NUM_THREADS %NUM_THREADS% -env PV_PLUGIN_PATH %PV_PLUGIN_PATH% -hosts 1 localhost -cores %MPI_PROCESSES% VestecCatalystEmulator.exe %DATA_PATH% .vtk 0.02 scripts\criticalPoints_Evaluation.py
@echo on
