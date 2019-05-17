#include "VestecSeedingAlgorithm.h"

#include <vtkPoints.h>
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
#include <random>
#include <cmath>

vtkStandardNewMacro(VestecSeedingAlgorithm);

//----------------------------------------------------------------------------
VestecSeedingAlgorithm::VestecSeedingAlgorithm()
{
  this->SetNumberOfInputPorts( 1 );
  this->SetNumberOfOutputPorts( 1 );
}

//----------------------------------------------------------------------------
VestecSeedingAlgorithm::~VestecSeedingAlgorithm()
{
}

//----------------------------------------------------------------------------
void VestecSeedingAlgorithm::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecSeedingAlgorithm::GetOutput()
{
  return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecSeedingAlgorithm::GetOutput(int port)
{
  return dynamic_cast<vtkDataSet*>(this->GetOutputDataObject(port));
}

//----------------------------------------------------------------------------
void VestecSeedingAlgorithm::SetOutput(vtkDataObject* d)
{
  this->GetExecutive()->SetOutputData(0, d);
}

//----------------------------------------------------------------------------
vtkDataObject* VestecSeedingAlgorithm::GetInput()
{
  return this->GetInput(0);
}

//----------------------------------------------------------------------------
vtkDataObject* VestecSeedingAlgorithm::GetInput(int port)
{
  return this->GetExecutive()->GetInputData(port, 0);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecSeedingAlgorithm::GetLabelHierarchyInput(int port)
{
  return dynamic_cast<vtkDataSet*>(this->GetInput(port));
}

//----------------------------------------------------------------------------
int VestecSeedingAlgorithm::ProcessRequest(vtkInformation* request,
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
int VestecSeedingAlgorithm::FillOutputPortInformation(
    int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int VestecSeedingAlgorithm::FillInputPortInformation(
                                               int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}


int VestecSeedingAlgorithm::RequestDataObject(vtkInformation* vtkNotUsed(request),
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
int VestecSeedingAlgorithm::RequestInformation(
                                         vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
                                      vtkInformationVector* vtkNotUsed(outputVector))
{
  // do nothing let subclasses handle it
  return 1;
}

//----------------------------------------------------------------------------
int VestecSeedingAlgorithm::RequestUpdateExtent(
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
int VestecSeedingAlgorithm::RequestData(
                                  vtkInformation* vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector* outputVector )
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkPolyData* output = dynamic_cast<vtkPolyData*>(outInfo->Get( vtkDataObject::DATA_OBJECT() ) );

  vtkInformation *inInfoGrid = inputVector[0]->GetInformationObject(0);
  vtkDataSet *inputGrid = dynamic_cast<vtkDataSet*>(inInfoGrid->Get(vtkDataObject::DATA_OBJECT()));


  int numPoints = inputGrid->GetNumberOfPoints();

  //Resulting points
  vtkPoints* pPoints = vtkPoints::New();
  pPoints->Allocate(numPoints * NumberOfPointsAroundSeed);

  // Create a c++11 random number generator
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

  for(int p = 0; p < numPoints; ++p)
  {
    //Get the initial position of the point
    double pos[3];
    inputGrid->GetPoint(p, pos);

    //Get lengt of each domain axis
    double bounds[6];
    inputGrid->GetBounds(bounds);

    //Calculate offset in each dimension
    double l_x = std::fabs((bounds[1] - bounds[0])) * (PercentOfDomain / 100);
    double l_y = std::fabs((bounds[3] - bounds[2])) * (PercentOfDomain / 100);
    double l_z = std::fabs((bounds[5] - bounds[4])) * (PercentOfDomain / 100);

    //Calculate distribution range for each axis
    std::uniform_real_distribution<> dis_x(pos[0] - l_x, pos[0] + l_x);
    std::uniform_real_distribution<> dis_y(pos[1] - l_y, pos[1] + l_y);
    std::uniform_real_distribution<> dis_z(pos[2] - l_z, pos[2] + l_z);

    for(int n = 0; n < NumberOfPointsAroundSeed; ++n)
    {
      double new_x = dis_x(gen);
      double new_y = dis_y(gen);
      double new_z = dis_z(gen);

      pPoints->InsertNextPoint(new_x, new_y, new_z);
    }
  }
  
  output->SetPoints(pPoints);
  return 1;
}

//----------------------------------------------------------------------------
void VestecSeedingAlgorithm::SetInput(vtkDataObject* input)
{
  this->SetInput(0, input);
}

//----------------------------------------------------------------------------
void VestecSeedingAlgorithm::SetInput(int index, vtkDataObject* input)
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
void VestecSeedingAlgorithm::AddInput(vtkDataObject* input)
{
  this->AddInput(0, input);
}

//----------------------------------------------------------------------------
void VestecSeedingAlgorithm::AddInput(int index, vtkDataObject* input)
{
  if(input)
  {
    this->AddInputDataObject(index, input);
  }
}