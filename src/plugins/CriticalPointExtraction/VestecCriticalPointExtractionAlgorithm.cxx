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

  //Set active array
  vtkDataArray *inScalars = this->GetInputArrayToProcess(0, inputVector);
  input->GetPointData()->SetActiveVectors(inScalars->GetName());

  //Compute critical pints
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
  std::cout << "Checking " << cells_num << " cells for critical points " << std::endl;

  unsigned long cp = 0;
  //Check for every cell if a critical point exists
  for(vtkIdType i=0; i < cells_num; i++) {
    if(PointInCell(input->GetCell(i), input)) {
        //add cell i to output
		cp++;
    }
  }
  std::cout << "Critical points found: " << cp << std::endl;
}


//int basic_isort3(int *a, int *b, int *c)
///* Input/Output: a, b, c. */
///* Sorts (&a,&b,&c) with <= 3 comparisons and <= 3 swaps. */
//{
//	int swaps = 0, aux;
//	if (*a > *b)
//		swap(*a, *b);
//	if (*b > *c)
//	{
//		swap(*b, *c);
//		if (*a > *b)
//			swap(*a, *b);
//	}
//	return (swaps);
//}

int CriticalPointExtractor::Sort3(vtkSmartPointer<vtkIdList> ids) 
{
	unsigned int tmp;
	unsigned int swaps = 0;
	if (ids->GetId(0) > ids->GetId(1))
	{
		tmp = ids->GetId(0);
		ids->SetId(0, ids->GetId(1));
		ids->SetId(1, tmp);
		swaps++;
	}

	if (ids->GetId(1) > ids->GetId(2))
	{
		tmp = ids->GetId(1);
		ids->SetId(1, ids->GetId(2));
		ids->SetId(2, tmp);
		swaps++;

		if (ids->GetId(0) > ids->GetId(1))
		{
			tmp = ids->GetId(0);
			ids->SetId(0, ids->GetId(1));
			ids->SetId(1, tmp);
			swaps++;
		}
	}
	return swaps;
}

int CriticalPointExtractor::Sort2(vtkSmartPointer<vtkIdList> ids)
{
	unsigned int tmp;
	unsigned int swaps = 0;
	if (ids->GetId(0) > ids->GetId(1))
	{
		tmp = ids->GetId(0);
		ids->SetId(0, ids->GetId(1));
		ids->SetId(1, tmp);
		swaps++;
	}
	return swaps;
}

/// can we pass to Positive directly the determinant matrix instead of the cell?
bool CriticalPointExtractor::Positive(vtkSmartPointer<vtkIdList> ids, vtkSmartPointer<vtkDataArray> vectors, unsigned long pertubationID)
{
	bool positivDir = false;

	// 1. Sort and check swap operations (check)
    // TODO: More generic version required. How to handle per pertubation
	int swapOperations = 0;
	if (ids->GetNumberOfIds() == 3) //TRIANGLES
		swapOperations = Sort3(ids);

	//create an eigen matrix
	MatrixXl vecMatrix(3, 3);
	for (vtkIdType i = 0; i < ids->GetNumberOfIds(); i++) {
		double * vecValues = vectors->GetTuple(ids->GetId(i));
		
		for (vtkIdType tuple = 0; tuple < vectors->GetNumberOfComponents(); tuple++) {
			//double to fixed precision
			vecMatrix(i, tuple) = toFixed(vecValues[tuple]);
		}
		vecMatrix(i,2) = 1;
	}

	//TODO: HACK 
	if (pertubationID != -1)
	{
		for (vtkIdType tuple = 0; tuple < vectors->GetNumberOfComponents() - 1; tuple++)
			vecMatrix(pertubationID, tuple) = 0;
	}
		
	//std::cout << " ######################################################################### " << std::endl;
	//std::cout << " " << vecMatrix(0, 0) << " " << vecMatrix(0, 1) << " " << vecMatrix(0, 2) << std::endl;
	//std::cout << " " << vecMatrix(1, 0) << " " << vecMatrix(1, 1) << " " << vecMatrix(1, 2) << std::endl;
	//std::cout << " " << vecMatrix(2, 0) << " " << vecMatrix(2, 1) << " " << vecMatrix(2, 2) << std::endl;
	//std::cout << " ######################################################################### " << std::endl;

	// 2. compute determinant sign
	unsigned long det = vecMatrix.determinant();
	if (det > 0) positivDir = true; else positivDir = false;

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		//positivDir *= -1;
		positivDir = false;
	}

	return positivDir;
}

bool CriticalPointExtractor::PointInCell(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid) {
	vtkSmartPointer<vtkIdList> ids = cell->GetPointIds();
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();

	// 1. compute the sign of the determinant of the cell
	// Get the sign (direction)
	bool initialDirection = Positive(ids, vectors);
	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0)
	for (int i = 0; i < ids->GetNumberOfIds(); i++) {
		// 2.2. compute the determinant sign again 
		// 2.3. check if it changes --> if so return false
		bool dir = Positive(ids, vectors, i);

		if (initialDirection != dir)
			return false;
	}
	return true; // the cell is critical, since the sign never change
}

unsigned long CriticalPointExtractor::toFixed(double val)
{
	//TODO: Some magic here
	val *= 100000000;
	return val * 100000000;
}