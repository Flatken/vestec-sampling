#ifndef VESTECADAPTOR_HEADER
#define VESTECADAPTOR_HEADER

#include <string>
#include <vector>

/**
 * This is a simple ParaView Catalyst adaptor which loads a time dependent dataset from disc and forwards it to the catalyst pipeline. This adaptor works also work in a
 * MPI application. When started in distributed (MPI) mode, it splits the datset according the number of preocesses. Therefore, each process only gets a portion of the total
 * domain. As input it gets a vector of filenames defining the simulation data to be read and a delta time describing the time between sthe snaphots in milliseconds
 */

#ifdef __cplusplus
extern "C"
{
#endif
	/**
	 * Create a Paraview Catalyst C/C++ pipeline
	 */
	void CatalystInitialize(std::vector<std::string> input, double deltaT, int numScripts, char* scripts[]);

	/**
	 * Clean shutdown of catalyst
	 */
	void CatalystFinalize();

	/**
	 * Load and split the vtkDataSet (loaded from disc) and forward that to our catalyst pipeline
	 */
	void CatalystCoProcess(double time, unsigned int timeStep, int lastTimeStep);
#ifdef __cplusplus
}
#endif
#endif