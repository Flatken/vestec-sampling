#include "VestecSamplingAlgorithm.h"

#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkCommand.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiProcessController.h>
#include <vtkDataSetWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkGeometryFilter.h>
#include <vtkDuplicatePolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <sstream>
#include <stdlib.h>

vtkStandardNewMacro(VestecSamplingAlgorithm);
vtkStandardNewMacro(UnsteadySource);

//----------------------------------------------------------------------------
VestecSamplingAlgorithm::VestecSamplingAlgorithm()
{
  this->SetNumberOfInputPorts( 2 );
  this->SetNumberOfOutputPorts( 1 );
  m_pUnsteadyData = UnsteadySource::New();

  double deltaT = 0.02;
	double t = 0;
  
	for(int i = 0; i < 500; i++)
	{
		timesteps.push_back(t);
		t += deltaT; 
	}
	m_pUnsteadyData->SetNumberOfTimeSteps(timesteps);

  m_pStreaklineFilter = vtkPStreaklineFilter::New();
  m_pStreaklineFilter->SetStaticMesh(1);
  m_pStreaklineFilter->SetStaticSeeds(0);
  m_pStreaklineFilter->SetForceReinjectionEveryNSteps(1);
  m_pStreaklineFilter->SetDisableResetCache(1);
  m_pStreaklineFilter->SetComputeVorticity(0);

  //vtkNew<vtkDataSetReader> pReader;
	//pReader->SetFileName("/unsecured/flat_ma/vestec/datasets/seeds/seeds_karman.vtk");
	//pReader->ReadAllScalarsOn();
	//pReader->Update();
  //m_pStreaklineFilter->SetInputData(1, pReader->GetOutput()); //Seeds
}

//----------------------------------------------------------------------------
VestecSamplingAlgorithm::~VestecSamplingAlgorithm()
{
  m_pStreaklineFilter->Delete();
}

//----------------------------------------------------------------------------
void VestecSamplingAlgorithm::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecSamplingAlgorithm::GetOutput()
{
  return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecSamplingAlgorithm::GetOutput(int port)
{
  return dynamic_cast<vtkDataSet*>(this->GetOutputDataObject(port));
}

//----------------------------------------------------------------------------
void VestecSamplingAlgorithm::SetOutput(vtkDataObject* d)
{
  this->GetExecutive()->SetOutputData(0, d);
}

//----------------------------------------------------------------------------
vtkDataObject* VestecSamplingAlgorithm::GetInput()
{
  return this->GetInput(0);
}

//----------------------------------------------------------------------------
vtkDataObject* VestecSamplingAlgorithm::GetInput(int port)
{
  return this->GetExecutive()->GetInputData(port, 0);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecSamplingAlgorithm::GetLabelHierarchyInput(int port)
{
  return dynamic_cast<vtkDataSet*>(this->GetInput(port));
}

//----------------------------------------------------------------------------
int VestecSamplingAlgorithm::ProcessRequest(vtkInformation* request,
                                     vtkInformationVector** inputVector,
                                     vtkInformationVector* outputVector)
{
  // Create an output object of the correct type.
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
  {
    return this->RequestDataObject(request, inputVector, outputVector);
  }
  // generate the data
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
  {
    return this->RequestData(request, inputVector, outputVector);
  }

  if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
  {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }

  // execute information
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
  {
    return this->RequestInformation(request, inputVector, outputVector);
  }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
int VestecSamplingAlgorithm::FillOutputPortInformation(
    int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int VestecSamplingAlgorithm::FillInputPortInformation(
                                               int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}


int VestecSamplingAlgorithm::RequestDataObject(vtkInformation* vtkNotUsed(request),
                                         vtkInformationVector** vtkNotUsed(inputVector),
         vtkInformationVector* outputVector )
    {
//RequestDataObject (RDO) is an earlier pipeline pass.
//During RDO, each filter is supposed to produce an empty data object of the proper type

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataSet* output = dynamic_cast<vtkDataSet*>(
    outInfo->Get( vtkDataObject::DATA_OBJECT() ) );

  if ( ! output )
  {
    output = vtkPolyData::New();
    outInfo->Set( vtkDataObject::DATA_OBJECT(), output );
    output->FastDelete();

    this->GetOutputPortInformation(0)->Set(
      vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType() );
  }

  return 1;
    }


//----------------------------------------------------------------------------
int VestecSamplingAlgorithm::RequestInformation(
                                         vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
                                      vtkInformationVector* vtkNotUsed(outputVector))
{
  // do nothing let subclasses handle it
  return 1;
}

//----------------------------------------------------------------------------
int VestecSamplingAlgorithm::RequestUpdateExtent(
                                          vtkInformation* vtkNotUsed(request),
    vtkInformationVector** inputVector,
    vtkInformationVector* vtkNotUsed(outputVector))
{
  int numInputPorts = this->GetNumberOfInputPorts();
  for (int i=0; i<numInputPorts; i++)
  {
    int numInputConnections = this->GetNumberOfInputConnections(i);
    for (int j=0; j<numInputConnections; j++)
    {
      vtkInformation* inputInfo = inputVector[i]->GetInformationObject(j);
      inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);
    }
  }
  return 1;
}

//----------------------------------------------------------------------------
// This is the superclasses style of Execute method.  Convert it into
// an imaging style Execute method.
int VestecSamplingAlgorithm::RequestData(
                                  vtkInformation* vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector* outputVector )
{
  vtkMultiProcessController *ctrl = vtkMultiProcessController::GetGlobalController();
	int mpiRank  	= ctrl->GetLocalProcessId();
  int numProcs  = ctrl->GetNumberOfProcesses();
  //Later on RequestData (RD) happens.
  //During RD each filter examines any inputs it has, then fills in that empty data object with real data.

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkDataSet* output = dynamic_cast<vtkDataSet*>(outInfo->Get( vtkDataObject::DATA_OBJECT() ) );

  vtkInformation *inInfoGrid = inputVector[0]->GetInformationObject(0);
  vtkDataSet *inputGrid = dynamic_cast<vtkDataSet*>(inInfoGrid->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *inInfoSeeds = inputVector[1]->GetInformationObject(0);
  vtkDataSet *inputSeeds = dynamic_cast<vtkDataSet*>(inInfoSeeds->Get(vtkDataObject::DATA_OBJECT()));


  ////TODO: Collect input seeds from all processes
  //if(mpiRank != 0)
  //{
  //  ctrl->Send(inputSeeds, 0, 999);
  //}else{
  //  for(int p = 0; p < numProcs; ++p)
  //  {
  //    vtmpiRankkPolyData* pd = vtkPolyData::New();
  //    controller->Receive(pd, i, 999);
  //    app->AddInputData(pd);
  //    pd->Delete();
  //  }
  //}
  
  double* d = inInfoGrid->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  int n  = inInfoGrid->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  double t = inInfoGrid->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  //set the active vector
  inputGrid->GetPointData()->SetActiveVectors("Vec");

  //Set timings for the unsteaday source
  timesteps.clear();
  timesteps.resize(n);
  inInfoGrid->Get( vtkStreamingDemandDrivenPipeline::TIME_STEPS(),&timesteps[0] );

  m_pUnsteadyData->InsertDataObject(t, inputGrid);

  std::vector<int> vecExtents(6);
  inInfoGrid->Get( vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),&vecExtents[0] );
  m_pUnsteadyData->SetExtent( vecExtents[0],
                              vecExtents[1],
                              vecExtents[2],
                              vecExtents[3],
                              vecExtents[4],
                              vecExtents[5]);
  
  if(timestepindex < 1)
  {
    timestepindex++;
    return 1;
  }

  //Collect seeds from every process since the tracer needs a consistent number of seeds accross MPI processes
  vtkUnstructuredGrid * pInputSeeds = vtkUnstructuredGrid::SafeDownCast(inputSeeds);
  vtkPoints * pSend = vtkPoints::New();
  vtkPoints * pRecv = vtkPoints::New();

  if(pInputSeeds->GetNumberOfPoints() > 0)
  {
    pSend->ShallowCopy(pInputSeeds->GetPoints());
  }
  ctrl->AllGatherV(pSend->GetData(), pRecv->GetData());

  vtkUnstructuredGrid * pCollectedSeeds = vtkUnstructuredGrid::New();
  pCollectedSeeds->SetPoints(pRecv);
  m_pStreaklineFilter->SetInputData(1, pCollectedSeeds); //Seeds
  m_pStreaklineFilter->SetInputConnection(0,m_pUnsteadyData->GetOutputPort()); //Grid 
  m_pStreaklineFilter->SetStartTime(t - 0.02);
  m_pStreaklineFilter->SetTerminationTime(t);
  m_pStreaklineFilter->Update();
  std::cout << "############## Integrating done for t: " << (t - 0.02) << " to " << t << std::endl; 
  std::cout << "###################################################################" << std::endl; 
  
  output->ShallowCopy(m_pStreaklineFilter->GetOutput());
  timestepindex++;
  return 1;
}

//----------------------------------------------------------------------------
void VestecSamplingAlgorithm::SetInput(vtkDataObject* input)
{
  this->SetInput(0, input);
}

//----------------------------------------------------------------------------
void VestecSamplingAlgorithm::SetInput(int index, vtkDataObject* input)
{
  if(input)
  {
    this->SetInputDataObject(index, input);
  }
  else
  {
    // Setting a NULL input removes the connection.
    this->SetInputDataObject(index, 0);
  }
}

//----------------------------------------------------------------------------
void VestecSamplingAlgorithm::AddInput(vtkDataObject* input)
{
  this->AddInput(0, input);
}

//----------------------------------------------------------------------------
void VestecSamplingAlgorithm::AddInput(int index, vtkDataObject* input)
{
  if(input)
  {
    this->AddInputDataObject(index, input);
  }
}