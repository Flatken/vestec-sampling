@echo off
set PV_PLUGIN_PATH=.;paraview-5.8/plugins/VestecPlugins

start smpd -d 3
mpiexec -env PV_PLUGIN_PATH %PV_PLUGIN_PATH% -hosts 1 localhost -cores 1 VestecCatalystEmulator.exe "D:\\vr_data\\KarmanVortexStreet\\2D\\" .vtk 0.02 criticalPoints_Evaluation.py
@echo on
