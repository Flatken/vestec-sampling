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
int VestecCriticalPointExtractionAlgorithm::FillOutputPortInformation(
	int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
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
  // get the input and output
  vtkDataSet* input = vtkDataSet::GetData(inputVector[0], 0);
  vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

  output->Print(std::cout);
  //CriticalPointExtractor cp_extractor;
  //cp_extractor.identify_critical_points(inputGrid);
  /// to-do MPI implementation
  ///  

  return 1;
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
    if(PointInCell(grid->GetCell(i),grid)) {
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
  std::vector<double*> points;
  for (vtkIdType i=0; i<ids->GetNumberOfIds(); i++) {
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