#include "SamplingPipeline.hpp"
#include "vtkSmartPointer.h"
#include "vtkMPIController.h"
#include "vtkMultiProcessController.h"
#include "vtkCPInputDataDescription.h"
#include "vtkObjectFactory.h"
#include <vtkDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <sstream>

vtkStandardNewMacro(SamplingPipeline);
vtkStandardNewMacro(CPTimeSource);

SamplingPipeline::SamplingPipeline()
{
    this->m_pSeeds = nullptr;
	vtkNew<vtkDataSetReader> pReader;
	pReader->SetFileName("/unsecured/flat_ma/vestec/datasets/seeds/seeds_karman.vtk");
	pReader->ReadAllScalarsOn();
	pReader->Update();

	this->m_pSeeds = pReader->GetPolyDataOutput()->NewInstance()vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
	int mpiRank  	= ctrl->GetLocalProcessId();
	this->m_pSeeds->ShallowCopy(pReader->GetPolyDataOutput());

	this->m_pInput = CPTimeSource::New();
	this->m_pStreakLineFilter = vtkPStreaklineFilter::New();
    m_pStreakLineFilter->SetStaticMesh(1);
    m_pStreakLineFilter->SetStaticSeeds(0);
    m_pStreakLineFilter->SetInputConnection(0,this->m_pInput->GevtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
	int mpiRank  	= ctrl->GetLocalProcessId();
	m_pStreakLineFilter->SetInputData(1,this->m_pSeeds); //SeedsvtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
	int mpiRank  	= ctrl->GetLocalProcessId();
	m_pStreakLineFilter->SetForceReinjectionEveryNSteps(1);
}



SamplingPipeline::~SamplingPipeline()
{
	std::cout << "[SamplingPipeline::~SamplingPipeline()]" << std::endl;
}

void SamplingPipeline::InitTimingInformation(int numTimeSteps, double deltaT)
{
	this->deltaT = deltaT;
	this->timesteps.clear();  
	double t = 0;
	for(int i = 0; i < numTimeSteps; i++)
	{
		this->timesteps.push_back(t);
		t += deltaT; 
	}
	this->m_pInput->SetNumberOfTimeSteps(this->timesteps);
	
}

int SamplingPipeline::RequestDataDescription(vtkCPDataDescription *dataDescription)
{
	if (!dataDescription)
	{
		vtkWarningMacro("data description is NULL");
		return 0;
	}
	return 1;
}

int SamplingPipeline::CoProcess(vtkCPDataDescription *dataDescription)
{
	vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
	int mpiRank  	= ctrl->GetLocalProcessId();
	int timestep 	= dataDescription->GetTimeStep();
	double dTime  	= dataDescription->GetTime();

	std::stringstream str;
    str << "Filename_rank_" << mpiRank << "_"<< timestep << ".vtk";

	//Register new timestep to our time dependent source
	vtkImageData *outImage = vtkImageData::SafeDownCast(dataDescription->GetInputDescriptionByName("input")->GetGrid());
    
	//Activate correct vector
	outImage->GetPointData()->SetActiveVectors("Vec");
	
	this->m_pInput->InsertDataObject(dTime, dataDescription->GetInputDescriptionByName("input")->GetGrid());
	
	//vtkImageData *outImage = vtkImageData::SafeDownCast(dataDescription->GetInputDescriptionByName("input")->GetGrid());

	if(timestep==0) //We need two timesteps available in input source
		return 1;

	//Read correct seeds
	//Todo: Later passed by catalyst pipeline
	std::stringstream seedName;
    seedName << "/unsecured/flat_ma/vestec/datasets/seeds/seeds_ttk_critical_points_vorticity_"<< timestep << ".vtk";
	vtkNew<vtkDataSetReader> pReader;
	pReader->SetFileName(seedName.str().c_str());
	pReader->ReadAllScalarsOn();
	pReader->Update();
	
	
	if(pReader->GetPolyDataOutput()->GetNumberOfPoints() > 0)
	{
		if(this->m_pSeeds)
			this->m_pSeeds->Delete();

		this->m_pSeeds = pReader->GetPolyDataOutput()->NewInstance();
		this->m_pSeeds->ShallowCopy(pReader->GetPolyDataOutput());
		//this->m_pSeeds->Print(std::cout);
		m_pStreakLineFilter->SetInputData(1,this->m_pSeeds); //Seeds
	}
	
	m_pStreakLineFilter->SetStartTime(dTime - this->deltaT);
	m_pStreakLineFilter->SetTerminationTime(dTime);
	//TODO: Check without that line the computation does not work
	m_pStreakLineFilter->SetDisableResetCache(1);
	m_pStreakLineFilter->Update();
	vtkPolyData* out = m_pStreakLineFilter->GetOutput();
    
	if(mpiRank == 0)
		std::cout << "[SamplingPipeline::CoProcess] Rank: " <<  mpiRank << " Time: " << dTime << " Timestep: " << timestep << " Points in output: " << out->GetNumberOfPoints()<< std::endl;
    
	//Write streaklines to disc
	vtkDataSetWriter* pWriter = vtkDataSetWriter::New();
	pWriter->SetInputData(out);
	pWriter->SetFileTypeToBinary();
    pWriter->SetFileName(str.str().c_str());
	pWriter->Update();
	pWriter->Delete();

	//this->m_pInput->Delete();
	
	return 1;
}

int SamplingPipeline::Finalize()
{
	std::cout << "[SamplingPipeline::Finalize()]" << std::endl;
	return 1;
}