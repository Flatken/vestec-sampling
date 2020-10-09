@echo off
set PATH=.;paraview-5.8/plugins/VestecPlugins;plugins/TopologyToolKit;%PATH%
set PYTHONPATH=Lib/site-packages
set PV_PLUGIN_PATH=.;paraview-5.8/plugins/VestecPlugins;plugins/TopologyToolKit;

start smpd -d 3
mpiexec -env OMP_NUM_THREADS %2 -env OMP_DISPLAY_ENV true -env OMP_PROC_BIND spread -env OMP_PLACES cores -hosts 1 localhost -cores %1 VestecCatalystEmulator.exe "C:\VESTEC\Datasets\lorenz" .vtk 1 criticalPoints_Lorenz.py
@echo on
