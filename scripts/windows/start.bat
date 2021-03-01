@echo off
set PV_PLUGIN_PATH=.;paraview-5.8/plugins/VestecPlugins

set OMP_NUM_THREADS=1
start smpd -d 3
mpiexec -env OMP_NUM_THREADS 1 -env PV_PLUGIN_PATH %PV_PLUGIN_PATH% -hosts 1 localhost -cores 1 VestecCatalystEmulator.exe "D:\\vr_data\\KarmanVortexStreet\\2D\\" .vtk 0.02 scripts\criticalPoints_Evaluation.py
@echo on
