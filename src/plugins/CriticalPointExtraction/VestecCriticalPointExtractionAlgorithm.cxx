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

  CriticalPointExtractor cp_extractor;
  cp_extractor.identify_critical_points(input, output);
  /// to-do MPI implementation
  ///  

  return 1;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::identify_critical_points(	vtkSmartPointer<vtkDataSet> input,
																				vtkSmartPointer<vtkDataSet> output) {
  vtkIdType cells_num = input->GetNumberOfCells();
 
  //Check for every cell if a critical point exists
  for(vtkIdType i=0; i<cells_num; i++) {
	
    if(PointInCell(input->GetCell(i), input)) {
      // add cell i to output
    }
  }
}

/// can we pass to Positive directly the determinant matrix instead of the cell?
bool CriticalPointExtractor::Positive(MatrixXl matrix)
{
	bool positivDir = false;
	// 1. TODO sort (check)
	// 2. compute determinant sign
	unsigned long det = matrix.determinant();
	if (det > 0) positivDir = true; else positivDir = false;
	// 3. check the number of swap operation while sorting


	return positivDir;
}

bool CriticalPointExtractor::PointInCell(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid) {
	vtkSmartPointer<vtkIdList> ids = cell->GetPointIds();
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();
	// 1. compute the sign of the determinant of the cell
	// 1.1 Create the eigen vector matrix
	MatrixXl vecMatrix(4, ids->GetNumberOfIds());
 
	for (vtkIdType i = 0; i < ids->GetNumberOfIds(); i++) {
		double* vecValues = vectors->GetTuple(i);
		for (vtkIdType tuple = 0; tuple < vectors->GetNumberOfComponents(); tuple++) {
			//double to fixed precision
			vecMatrix << toFixed(vecValues[tuple]);
		}
		vecMatrix << 1;
	}
	// Get the sign (direction)
	bool initialDirection = Positive(vecMatrix);

	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0)
	for (int i = 0; i < ids->GetNumberOfIds(); i++) {
		MatrixXl pertubationMatrix = vecMatrix;
		for (int tuple = 0; tuple < vectors->GetNumberOfComponents(); tuple++) 
		{
			pertubationMatrix(i, tuple) = 0;
		}

		// 2.2. compute the determinant sign again 
		// 2.3. check if it changes --> if so return false
		bool dir = Positive(pertubationMatrix);

		if (initialDirection != dir)
			return false;
	}
	return true; // the cell is critical, since the sign never change
}

unsigned long CriticalPointExtractor::toFixed(double val)
{
	char str[50];
	sprintf(str, "%0.*f", 24, val);
	return atof(str);
}