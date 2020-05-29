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
#include <vtkPointData.h>

#include <sstream>
#include <cmath>

#define M_PI 3.14159

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

  int localNumberOfSeeds = NumberOfPointsAroundSeed;
  int numPoints = inputGrid->GetNumberOfPoints();

  //Resulting points
  vtkPoints* pPoints = vtkPoints::New();
  pPoints->Allocate(numPoints * NumberOfPointsAroundSeed);

  // Create a c++11 random number generator
  std::random_device rd;
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

  //For every input point generate n seeds
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

    double l = 0;
    if(l_x != 0 && l_y != 0 && l_z != 0)
      l = std::min(std::min(l_x, l_y), l_z);
    else
      l = std::min(l_x, l_y);
    
    AbstractGenerator* genAngleY = nullptr;
    AbstractGenerator* genAngleZ = nullptr;
    AbstractGenerator* genLenght = nullptr;

    if(DistributionMode == 0)
    {
      genAngleY = new RandomNumberGenerator<std::uniform_real_distribution<>>(0,M_PI*2);
      genAngleZ = new RandomNumberGenerator<std::uniform_real_distribution<>>(0,M_PI*2);
      genLenght = new RandomNumberGenerator<std::uniform_real_distribution<>>(0,l);
    }else if(DistributionMode == 1)
    {
      genAngleY = new RandomNumberGenerator<std::normal_distribution<>>(0,M_PI*2);
      genAngleZ = new RandomNumberGenerator<std::normal_distribution<>>(0,M_PI*2);
      genLenght = new RandomNumberGenerator<std::normal_distribution<>>(0,l);
    }else{
      std::cout << "Error: Unknown distribution mode" << std::endl;
    }

    //Detetmine the number of seeds with the active scalar value and their range
    if(UseScalarRange)
    {
      vtkDataArray *inScalars = this->GetInputArrayToProcess(0, inputVector);
      double range[2];
      inScalars->GetRange(range);
      std::cout << "Name " << inScalars->GetName() << " Range "<<  range[0] << "  " << range[1] << std::endl;
      localNumberOfSeeds = 0;
    }

    //Generate the random seeds and add them to the output point array
    for(int n = 0; n < localNumberOfSeeds; ++n)
    {
      double distance = genLenght->gen();
      double angleY   = genAngleY->gen();
      double angleZ   = genAngleZ->gen();

      double new_pos[3];
      if(l_x != 0 && l_y != 0 && l_z != 0)
      {
        new_pos[0] = pos[0] + distance * std::cos(angleZ) * std::sin(angleY);
        new_pos[1] = pos[1] + distance * std::sin(angleZ);
        new_pos[2] = pos[2] + distance * std::cos(angleZ) * std::cos(angleY);
      }else{
        new_pos[0] = pos[0] + distance * std::cos(angleY);
        new_pos[1] = pos[1] + distance * std::sin(angleY);
        new_pos[2] = pos[2];
      }
      pPoints->InsertNextPoint(new_pos[0], new_pos[1], new_pos[2]);
    }

    delete genLenght;
    delete genAngleY;
    delete genAngleZ;
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