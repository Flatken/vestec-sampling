## Description
Evaluation tool for the execution of a Python based catalyst pipeline. The tool emulates a running simulation by loading a time dpendent dataset timestep for timestep from disc and forwards these data to the paraview catalyst pipeline for execution. The tool supports the evaluation of filters using distributed memory parallelizm by splitting each dataset to the given number of MPI processes. Each MPI process then gets only a portion (chunk) of the total domain. Furthermore, the repository includes own Paraview plugins with VTK filters for data sampling of the raw data. 

The idea is to use topological information extracted using the Topology Toolkit (TTK) together with a pathline snippet tracing algorithm to sample the original raw data. From these information it should be possible to reconstruct an approximated dataset.

Requirements to build the Software:

- GCC (Linux) or MSVC (Windows) Compiler with C++17 support
- MPI installation of your choise (MVAPICH2, Microsoft MPI on Windows)
- CMake >= 3.12
- Python >= 3.4
- Boost >= 1.69


The tool has the following external dependencies integrated as git submodules. You do not have to install them separately:

- [TTK (0.9.7)](https://topology-tool-kit.github.io/)
- [ParaView (5.6)](https://www.paraview.org/)


The software is funded by the [VESTEC EU Project](https://www.vestec-project.eu/ "VESTEC EU Project")

## Build the software

Clone the git repository to your local hard drive. Choose a folder ($FOLDER)
```
cd $FOLDER
git clone https://github.com/Flatken/vestec-sampling.git src
```

Update and initialize the git submodules recursivly. 

```
cd src
git submodule update --init --recursive 
```
Then, patch Paraview in order make it compatible with TTK.
```
./patch_paraview.sh
```

Build the external dependencies by using the folowing shell script (ParaView and TTK)

```
cd $FOLDER
src/make_externals.sh
```

Build the evaluation tool itself uby sing the folowing shell script

```
cd $FOLDER
src/make.sh
```

## User Guide
- Exporting a ParaView Catalyst pipeline
	- See TKK tutorials on Catalyst 
- Requirements
	- At the moment only VTK\_STRUCTURED\_POINTS data in legacy format is supported as input data
	- The data naming for an unsteady dataset needs to be in the following format. Each file needs to be placed in the same directory. 
		- snapshot_??.vtk
			- Where ?? is the timestep with leading zeros
				- snapshot_00.vtk
				- snapshot_01.vtk
		
- Call the evaluation tool. You can also have look to the start.sh shell script
```
	cd $FOLDER/install/linux-release/bin
	mpiexec -np 1 ./VestecCatalystEmulator [DATA_PATH] [TIME STEPSIZE/DELTA] [CATALYST SCRIPT]
```
## Hints when using Phyton version 3.4
There seems to be a bug in paraview which leads to crashes when executing the application. You have to fix the following file
and then call src/make_release.sh again.

FIX: 
```
vim $FOLDER/install/linux-externals/lib/python3.4/site-packages/paraview/vtk.py
```
Change the following line by replacing version 3,4 with version 3,5!
```
if sys.version_info < (3,4):
```
## Support
If you have any questions or comments just write an e-mail to: markus.flatken@dlr.de
