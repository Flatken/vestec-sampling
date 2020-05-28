#include "VestecCriticalPointExtractionAlgorithm.h"

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

vtkStandardNewMacro(VestecCriticalPointExtractionAlgorithm);

//----------------------------------------------------------------------------
VestecCriticalPointExtractionAlgorithm::VestecCriticalPointExtractionAlgorithm()
{
  this->SetNumberOfInputPorts( 1 );
  this->SetNumberOfOutputPorts( 1 );
}

//----------------------------------------------------------------------------
VestecCriticalPointExtractionAlgorithm::~VestecCriticalPointExtractionAlgorithm()
{
}

//----------------------------------------------------------------------------
void VestecCriticalPointExtractionAlgorithm::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecCriticalPointExtractionAlgorithm::GetOutput()
{
  return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecCriticalPointExtractionAlgorithm::GetOutput(int port)
{
  return dynamic_cast<vtkDataSet*>(this->GetOutputDataObject(port));
}

//----------------------------------------------------------------------------
void VestecCriticalPointExtractionAlgorithm::SetOutput(vtkDataObject* d)
{
  this->GetExecutive()->SetOutputData(0, d);
}

//----------------------------------------------------------------------------
vtkDataObject* VestecCriticalPointExtractionAlgorithm::GetInput()
{
  return this->GetInput(0);
}

//----------------------------------------------------------------------------
vtkDataObject* VestecCriticalPointExtractionAlgorithm::GetInput(int port)
{
  return this->GetExecutive()->GetInputData(port, 0);
}

//----------------------------------------------------------------------------
vtkDataSet* VestecCriticalPointExtractionAlgorithm::GetLabelHierarchyInput(int port)
{
  return dynamic_cast<vtkDataSet*>(this->GetInput(port));
}

//----------------------------------------------------------------------------
int VestecCriticalPointExtractionAlgorithm::ProcessRequest(vtkInformation* request,
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
int VestecCriticalPointExtractionAlgorithm::FillOutputPortInformation(
    int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int VestecCriticalPointExtractionAlgorithm::FillInputPortInformation(
                                               int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}


int VestecCriticalPointExtractionAlgorithm::RequestDataObject(vtkInformation* vtkNotUsed(request),
                                         vtkInformationVector** vtkNotUsed(inputVector),
         vtkInformationVector* outputVector ) {
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
int VestecCriticalPointExtractionAlgorithm::RequestInformation(
                                         vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
                                      vtkInformationVector* vtkNotUsed(outputVector)) {
  // do nothing let subclasses handle it
  return 1;
}

//----------------------------------------------------------------------------
int VestecCriticalPointExtractionAlgorithm::RequestUpdateExtent(
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
int VestecCriticalPointExtractionAlgorithm::RequestData(
                                  vtkInformation* vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector* outputVector )
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkPolyData* output = dynamic_cast<vtkPolyData*>(outInfo->Get( vtkDataObject::DATA_OBJECT() ) );

  vtkInformation *inInfoGrid = inputVector[0]->GetInformationObject(0);
  vtkDataSet *inputGrid = dynamic_cast<vtkDataSet*>(inInfoGrid->Get(vtkDataObject::DATA_OBJECT()));

  identify_critical_points(inputGrid);
  /// to-do MPI implementation
  ///  

  // int localNumberOfSeeds = NumberOfPointsAroundSeed;
  // int numPoints = inputGrid->GetNumberOfPoints();

  // //Resulting points
  // vtkPoints* pPoints = vtkPoints::New();
  // pPoints->Allocate(numPoints * NumberOfPointsAroundSeed);
  
  // //For every input point generate n seeds
  // for(int p = 0; p < numPoints; ++p)
  // {
  //   //Get the initial position of the point
  //   double pos[3];
  //   inputGrid->GetPoint(p, pos);

  //   //Get lengt of each domain axis
  //   double bounds[6];
  //   inputGrid->GetBounds(bounds);

  //   //Calculate offset in each dimension
  //   double l_x = std::fabs((bounds[1] - bounds[0])) * (PercentOfDomain / 100);
  //   double l_y = std::fabs((bounds[3] - bounds[2])) * (PercentOfDomain / 100);
  //   double l_z = std::fabs((bounds[5] - bounds[4])) * (PercentOfDomain / 100);

  //   double l = 0;
  //   if(l_x != 0 && l_y != 0 && l_z != 0)
  //     l = std::min(std::min(l_x, l_y), l_z);
  //   else
  //     l = std::min(l_x, l_y);
    
  //   // AbstractGenerator* genAngleY = nullptr;
  //   // AbstractGenerator* genAngleZ = nullptr;
  //   // AbstractGenerator* genLenght = nullptr;

  //   if(DistributionMode == 0)
  //   {
  //     // genAngleY = new RandomNumberGenerator<std::uniform_real_distribution<>>(0,M_PI*2);
  //     // genAngleZ = new RandomNumberGenerator<std::uniform_real_distribution<>>(0,M_PI*2);
  //     // genLenght = new RandomNumberGenerator<std::uniform_real_distribution<>>(0,l);
  //   }else if(DistributionMode == 1)
  //   {
  //     // genAngleY = new RandomNumberGenerator<std::normal_distribution<>>(0,M_PI*2);
  //     // genAngleZ = new RandomNumberGenerator<std::normal_distribution<>>(0,M_PI*2);
  //     // genLenght = new RandomNumberGenerator<std::normal_distribution<>>(0,l);
  //   }else{
  //     std::cout << "Error: Unknown distribution mode" << std::endl;
  //   }

  //   //Detetmine the number of seeds with the active scalar value and their range
  //   if(UseScalarRange)
  //   {
  //     vtkDataArray *inScalars = this->GetInputArrayToProcess(0, inputVector);
  //     double range[2];
  //     inScalars->GetRange(range);
  //     std::cout << "Name " << inScalars->GetName() << " Range "<<  range[0] << "  " << range[1] << std::endl;
  //     localNumberOfSeeds = 0;
  //   }

  //   // //Generate the random seeds and add them to the output point array
  //   // for(int n = 0; n < localNumberOfSeeds; ++n)
  //   // {
  //   //   double distance = 0;//genLenght->gen();
  //   //   double angleY   = 0;//genAngleY->gen();
  //   //   double angleZ   = 0;//genAngleZ->gen();

  //   //   double new_pos[3];
  //   //   if(l_x != 0 && l_y != 0 && l_z != 0)
  //   //   {
  //   //     new_pos[0] = pos[0] + distance * std::cos(angleZ) * std::sin(angleY);
  //   //     new_pos[1] = pos[1] + distance * std::sin(angleZ);
  //   //     new_pos[2] = pos[2] + distance * std::cos(angleZ) * std::cos(angleY);
  //   //   }else{
  //   //     new_pos[0] = pos[0] + distance * std::cos(angleY);
  //   //     new_pos[1] = pos[1] + distance * std::sin(angleY);
  //   //     new_pos[2] = pos[2];
  //   //   }
  //   //   pPoints->InsertNextPoint(new_pos[0], new_pos[1], new_pos[2]);
  //   // }

  //   // delete genLenght;
  //   // delete genAngleY;
  //   // delete genAngleZ;
  // }
  // output->SetPoints(pPoints);
  // return 1;
}

//----------------------------------------------------------------------------
void VestecCriticalPointExtractionAlgorithm::SetInput(vtkDataObject* input)
{
  this->SetInput(0, input);
}

//----------------------------------------------------------------------------
void VestecCriticalPointExtractionAlgorithm::SetInput(int index, vtkDataObject* input)
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
void VestecCriticalPointExtractionAlgorithm::AddInput(vtkDataObject* input)
{
  this->AddInput(0, input);
}

//----------------------------------------------------------------------------
void VestecCriticalPointExtractionAlgorithm::AddInput(int index, vtkDataObject* input)
{
  if(input)
  {
    this->AddInputDataObject(index, input);
  }
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> CriticalPointExtractor::identify_critical_points(vtkSmartPointer<vtkDataSet> grid) {
  vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
  vtkIdType cells_num = grid->GetNumberOfCells();
  // bool is_critical = false;
  for(vtkIdType i=0; i<cells_num; i++) {
    // vtkSmartPointer<vtkTetra> t = vtkTetra::SafeDownCast(tet_mesh->GetCell(i));
    // float t_degree = compute_degree(t);
    // // here we need to understand how to discretize the critical points..
    // // likely we have to understand 
    if(PointInCell(grid->GetCell(i)),grid) {
      // add cell i to output
    }
  }
  return output;
}

/// can we pass to Positive directly the determinant matrix instead of the cell?
int CriticalPointExtractor::Positive(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid)
{
  vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();
  // 0. initialize vector array of cell points vectors
  vtkSmartPointer<vtkIdList> ids = cell->GetPointIds();
  std::vector<std::unique_ptr<double>> points;
  for (auto i : ids) {
    points.push_back(vectors->GetTuple(i));
  }
  // 1. --> convert to fixed precision (float to long)
  // 2. sort (check)
  // 3. create vector matrix and compute determinant sign
  // 4. check the number of swap operation while sorting

}

bool CriticalPointExtractor::PointInCell(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid) {
  // 1. compute the sign of the determinant of the cell  
  // 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
  // 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0)
  // 2.2. compute the determinant sign again 
  // 2.3. check if it changes --> if so return false
  return true; // the cell is critical, since the sign never change
}