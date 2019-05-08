## Description
Evaluation tool for the execution of a Python based catalyst pipeline. The tool emulates a running simulation by loading a time dpendent dataset timestep for timestep from disc and forwards these data to the paraview catalyst pipeline for execution. The tool supports the evaluation of filters using distributed memory parallelizm by splitting each dataset to the given number of MPI processes. Each MPI process then gets only a portion (chunk) of the total domain. Furthermore, the repository includes own Paraview plugins with VTK filters for data sampling of the raw data. 

The idea is to use topological information extracted using the Topology Toolkit (TTK) together with a pathline snippet tracing algorithm to sample the original raw data. From these information it should be possible to reconstruct an approximated dataset.

Prerequirements to build the Software:

- GCC (Linux) or MSVC (Windows) Compiler with C++17 support
- MPI installation of your choise (MVAPICH2, Microsoft MPI on Windows)
- CMake >= 3.12
- Python >= 3.4


The tool has the following external dependencies:

- [TTK (0.9.7)](https://topology-tool-kit.github.io/)
- [ParaView (5.6)](https://www.paraview.org/)

These dependecies are integrated as git submodules.

The software is funded by the [VESTEC EU Project](https://www.vestec-project.eu/ "VESTEC EU Project")

## Build the software

Clone the git repository to your local hard drive 
```
cd $FOLDER
git clone https://github.com/Flatken/vestec-sampling.git src
```

Update and initialize the git submodules recursivly. 

```
cd src
git submodule update --init --recursive 
```

Build the external dependencies (ParaView and TTK)

```
cd $FOLDER
src/make_externals.sh
```

Build the evaluation tool itself

```
cd $FOLDER
src/make_release.sh
```

## User Guide
- Exporting a ParaView Catalyst pipeline
	- See TKK tutorials on Catalyst 
- Requirements
	- At the moment only VTK\_STRUCTURED\_POINTS data in legacy format is supported
	- The data naming for an unsteady dataset needs to be in the following format. Each file needs to be placed in the same directory. 
		- snapshot_??.vtk
			- Where ?? is the timestep with leading zeros
				- snapshot_00.vtk
				- snapshot_01.vtk
		
- Calling the evaluation tool
	- ./start.sh [DATA_PATH] [TIME STEPSIZE/DELTA]
Following

## Support
If you have any questions or comments just write an e-mail to: markus.flatken@dlr.de
