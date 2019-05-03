#include <iostream>
#include "VestecAdaptor.hpp"
//#include "../CPFilters/SamplingPipeline.hpp"

#include <vtkType.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkDataObject.h>
#include <vtkDataSetReader.h>
#include <vtkUniformGridPartitioner.h>
#include <vtkMultiProcessController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkResampleWithDataSet.h>
#include <vtkCPPythonScriptPipeline.h>

namespace
{
	vtkCPProcessor* Processor = NULL;
	vtkDataObject* VTKGrid;
	std::vector<std::string> inputFiles;
}

void CatalystInitialize(std::vector<std::string> input, double deltaT, int numScripts, char* scripts[])
{
	//Store filenames for later usage
	inputFiles = input;

	if (Processor == NULL)
	{
		Processor = vtkCPProcessor::New();
		Processor->Initialize();
	}
	else
	{
		Processor->RemoveAllPipelines();
	}

	for(int i=3;i<numScripts;i++)
    {
    	vtkNew<vtkCPPythonScriptPipeline> pipeline;
    	pipeline->Initialize(scripts[i]);
    	Processor->AddPipeline(pipeline.GetPointer());
		std::cout << "Adding script " << scripts[i] << " to catalyst pipeline" << std::endl;
	}
	//Create our Catalyst pipeline
	//TODO: How to establish such a pipeline using the vtkCPPythonScriptPipeline with our filter? 
	//SamplingPipeline* pPathlineFilter = SamplingPipeline::New();
	//vtkCPPythonScriptPipeline* pPythonPipeline = vtkCPPythonScriptPipeline::New();
	//pPythonPipeline-
	//pPathlineFilter->InitTimingInformation(input.size(), deltaT);
	//Processor->AddPipeline(pPathlineFilter);
}

void CatalystFinalize()
{
	if (Processor)
	{
		Processor->Delete();
		Processor = NULL;
	}
	if (VTKGrid)
	{
		VTKGrid = NULL;
	}
}

void CatalystCoProcess( double time, unsigned int timeStep, int lastTimeStep)
{
	vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
	int mpiRank = ctrl->GetLocalProcessId();
	int mpiRanks = ctrl->GetNumberOfProcesses();

	vtkNew<vtkCPDataDescription> dataDescription;
	dataDescription->AddInput("input");

	//Set the timing information of the input data
	dataDescription->SetTimeData(time, timeStep);
	if (lastTimeStep == true)
	{
		// assume that we want to all the pipelines to execute if it
		// is the last time step.
		dataDescription->ForceOutputOn();
	}
	if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
	{
		//Load the correct data for the timeStep
		vtkDataSetReader* pReader = vtkDataSetReader::New();
		//Get the file for the correct timestep (files are here dataPath)
		pReader->SetFileName(inputFiles[timeStep].c_str());
		pReader->SetReadAllScalars(true);
		pReader->Update();

		//Check the current data type
		if (pReader->ReadOutputType() == VTK_STRUCTURED_POINTS)
		{
			//Split the domain and pass only a portion (block) to the catalyst pipeline
			vtkUniformGridPartitioner* pPartitioner = vtkUniformGridPartitioner::New();
			pPartitioner->SetInputData(pReader->GetOutput());
			pPartitioner->SetNumberOfPartitions(mpiRanks);
			pPartitioner->Update();
			
			//Pass the scalar fields to the portion (block) 
			//Todo: Check if there is a better way of passing scalars and vectors
			vtkMultiBlockDataSet *mbds = vtkMultiBlockDataSet::SafeDownCast(pPartitioner->GetOutput());
			vtkResampleWithDataSet* pResample = vtkResampleWithDataSet::New();
			pResample->SetInputData(mbds->GetBlock(mpiRank));
			pResample->SetSourceData(pReader->GetOutput());
			pResample->PassPointArraysOn();
			pResample->PassCellArraysOff();
			pResample->PassFieldArraysOff();
            pResample->Update();

			//Copy the portion for later processing
			VTKGrid = pResample->GetOutput()->NewInstance();
			VTKGrid->ShallowCopy(pResample->GetOutput());
			
			//Cleanup the filters used to split the dataset
			pPartitioner->Delete();
			pResample->Delete();
		}
		else if (pReader->ReadOutputType() == VTK_STRUCTURED_POINTS)
		{

		}
		else if (pReader->ReadOutputType() == VTK_STRUCTURED_POINTS)
		{

		}
		else {
			std::cout << "Unsupported vtk data format! " << std::endl;
		}
			
		//Pass grid (portion) to the catalyst pipeline
		dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
		Processor->CoProcess(dataDescription.GetPointer());
	

		//Cleanup memory
		pReader->Delete();
	}
}
