#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <algorithm>

#include "Adaptors/VestecAdaptor.hpp"


#include <experimental/filesystem>

using namespace std;
namespace filesys = std::experimental::filesystem;

/**
 * Returns a vector with the full path to each simulation file
 */
std::vector<std::string> getAllFilesInDir(const std::string &dirPath, const std::vector<std::string> fileSkipList = {})
{

	// Create a vector of string
	std::vector<std::string> listOfFiles;
	try {
		// Check if given path exists and points to a directory
		if (filesys::exists(dirPath) && filesys::is_directory(dirPath))
		{
			// Create a Recursive Directory Iterator object and points to the starting of directory
			filesys::recursive_directory_iterator iter(dirPath);

			// Create a Recursive Directory Iterator object pointing to end.
			filesys::recursive_directory_iterator end;

			// Iterate till end
			while (iter != end)
			{
				// Check if current entry is a directory and if exists in skip list
				if (filesys::is_directory(iter->path()))
				{
					// Skip the iteration of current directory pointed by iterator
#ifdef USING_BOOST
					// Boost Fileystsem  API to skip current directory iteration
					iter.no_push();
#else
					// c++17 Filesystem API to skip current directory iteration
					iter.disable_recursion_pending();
#endif
				}
				else
				{
					string::npos;
					if (std::find(fileSkipList.begin(), fileSkipList.end(), iter->path().extension()) != fileSkipList.end())
					{
						listOfFiles.push_back(iter->path().string());
					}	
				}

				error_code ec;
				// Increment the iterator to point to next entry in recursive iteration
				iter.increment(ec);
				if (ec) {
					std::cout << "Error While Accessing : " << iter->path().string() << " :: " << ec.message() << '\n';
				}
			}
		}
	}
	catch (std::system_error & e)
	{
		std::cout << "Exception :: " << e.what();
	}
	return listOfFiles;
}


int main(int argc, char** argv) {
	if (argc <= 3)
	{
		printf("Usage: VesteVestecSampling [PATH-TO-DATA] [DELTA_TIME] [CATALYST_PYTHON_SCRIPTS] \n");
		return 0;
	}
	
	// Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

	//get number of files (each file represents one timestep of the simulation)
	std::vector<std::string>files = getAllFilesInDir(std::string(argv[1]), { ".vtk" });
	std::sort(files.begin(), files.end());
	
	if(world_rank == 0)
	{
		std::cout << "######################################################################## " << std::endl;
		std::cout << "This an example programm for VESTEC to calculate proxies and do sampling " << std::endl;
    	std::cout << "Reading simulation files from: " << std::string(argv[1]) << std::endl;
		std::cout << "Number of input files: " << files.size() << std::endl;
		std::cout << "Delta between timesteps (ms): " << std::string(argv[2]) << std::endl;
		std::cout << "Catalyst Python script: " << std::string(argv[3]) << std::endl;
		std::cout << "######################################################################## " << std::endl;
	}

	unsigned int numberOfTimeSteps = files.size();
	unsigned int timeStep = 0;
	double time = 0;
	double stepsize = stod(std::string(argv[2]));

	#ifdef USE_CATALYST
		CatalystInitialize(files, stepsize, argc, argv);
	#endif

	for (timeStep = 0; timeStep < numberOfTimeSteps; timeStep++)
	{
		std::cout << "Execute pipeline for t: "<< time <<" and step: " << timeStep << std::endl;
#ifdef USE_CATALYST
		int lastTimeStep = 0;
		if (timeStep == numberOfTimeSteps - 1)
		{
			lastTimeStep = 1;
		}
	    CatalystCoProcess(time, timeStep, lastTimeStep);
		time += stepsize;
#endif
	}

#ifdef USE_CATALYST
	CatalystFinalize();
#endif

    // Finalize the MPI environment.
    MPI_Finalize();
}

