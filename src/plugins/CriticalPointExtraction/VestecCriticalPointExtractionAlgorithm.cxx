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
#include <vtkCellArray.h>
#include <vtkGenericCell.h>

#include <sstream>
#include <cmath>


#include <omp.h>

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

  //Compute critical points
  std::vector<double*> singularities;
  double p1[3] = {0, 0, 0};
  double p2[3] = {0.001, 0.001, 0};
  double p3[3] = {0.01, 0.01, 0};
  double p4[3] = {0.02, 0.02, 0};
  singularities.push_back(p1);
  singularities.push_back(p2);
  singularities.push_back(p3);
  singularities.push_back(p4);
  CriticalPointExtractor cp_extractor;
  cp_extractor.identify_critical_points(input, output, singularities);
  /// to-do MPI implementation
  ///  

  return 1;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::identify_critical_points(	vtkSmartPointer<vtkDataSet> input,
																				vtkSmartPointer<vtkDataSet> output, std::vector<double*> singlarities) {
	vtkSmartPointer<vtkPolyData> outputData = vtkPolyData::SafeDownCast(output);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  
	vtkIdType cells_num = input->GetNumberOfCells();
	std::cout << "Checking " << cells_num << " cells for critical points " << std::endl;

	//Set zero id to end 
	ZERO_ID = input->GetNumberOfPoints();

	//Set the singularity vector value (normally 0 in any dimension)
	currentSingularity[0] = singlarities[0][0];
	currentSingularity[1] = singlarities[0][1];
	currentSingularity[2] = singlarities[0][2];
#
	//Prepare for parallel computation
	int numThreads = 12;
	omp_set_num_threads(numThreads);
	std::vector<vtkSmartPointer<vtkGenericCell>> vecCellPerThread;
	for(int x=0; x < numThreads; x++)
		vecCellPerThread.push_back(vtkSmartPointer<vtkGenericCell>::New());
	
	//Check for every cell if a critical point exists
	#pragma omp parallel for 
	for(vtkIdType i = 0; i < cells_num; i++) {
		//get the cell for the current thread
		int threadIdx = omp_get_thread_num();
		input->GetCell(i, vecCellPerThread[threadIdx]);

		//If the cell contains a critical point add them to the output
		if(PointInCell(vecCellPerThread[threadIdx], input)) {
			vtkSmartPointer<vtkIdList> ids = vecCellPerThread[threadIdx]->GetPointIds();
			vtkSmartPointer<vtkIdList> new_ids = vtkSmartPointer<vtkIdList>::New();;
			for (int index = 0; index < ids->GetNumberOfIds(); index++)
			{
				double pCoords[3];
				vtkIdType newCellID;
				input->GetPoint(ids->GetId(index), pCoords);
#pragma omp critical
				newCellID = points->InsertNextPoint(pCoords[0], pCoords[1], pCoords[2]);
				new_ids->InsertNextId(newCellID);
			}
#pragma omp critical
			cells->InsertNextCell(new_ids);
		}
	}

  //Add points and cells to polydata
  outputData->SetPoints(points); 
  outputData->SetPolys(cells);
  std::cout << "Critical points found: " << cells->GetNumberOfCells() << std::endl;
}

bool CriticalPointExtractor::PointInCell(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid) {
	vtkSmartPointer<vtkIdList> ids = cell->GetPointIds();
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();

	//std::cout << " ################### Point in Cell ###################################################### " << std::endl;
	// 1. compute the sign of the determinant of the cell
	// Get the determinat and direction
	double initialDeterminant = Positive(ids, vectors);
	bool initialDirection     = DeterminatCounterClockWise(initialDeterminant);

	//Check for non data values (vector is zero) 
	if (initialDeterminant == 0)
		return false;

	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0)
	for (int i = 0; i < ids->GetNumberOfIds(); i++) {
		// 2.2. compute the determinant sign again 
		// 2.3. check if it changes --> if so return false
		double tmpDeterminat = Positive(ids, vectors, i);
		bool tmpDirection    = DeterminatCounterClockWise(tmpDeterminat);

		if (initialDirection != tmpDirection)
			return false;
	}
	//std::cout << " \t\t Cell is critical " << std::endl;
	return true; // the cell is critical, since the sign never change
}

/// can we pass to Positive directly the determinant matrix instead of the cell?
double CriticalPointExtractor::Positive(vtkSmartPointer<vtkIdList> ids, vtkSmartPointer<vtkDataArray> vectors, long pertubationID){
	vtkSmartPointer<vtkIdList> tmpIds = vtkSmartPointer<vtkIdList>::New();
	tmpIds->DeepCopy(ids);

	if (pertubationID != -1)
	{
		tmpIds->SetId(pertubationID,ZERO_ID);
	}

	// 1. Sort and check swap operations (check)
    // TODO: More generic version required. How to handle per pertubation
	int swapOperations = 0;
	if (tmpIds->GetNumberOfIds() == 3) //TRIANGLES(2D)
		swapOperations = Sort3(tmpIds);
	
	//create an eigen matrix
	MatrixXl vecMatrix;
	for (vtkIdType tuple = 0; tuple < tmpIds->GetNumberOfIds(); tuple++) {
		double vecValues[3];
		if (tmpIds->GetId(tuple) != ZERO_ID)
		{
			vectors->GetTuple(tmpIds->GetId(tuple), vecValues);
		}
		else
		{
			vecValues[0] = currentSingularity[0];
			vecValues[1] = currentSingularity[1];
			vecValues[2] = currentSingularity[2];
		}
		
		for (vtkIdType i = 0; i < vectors->GetNumberOfComponents(); i++) {
			//TODO: long long to fixed precision
			//vecMatrix(tuple, i) = vecValues[i];
			vecMatrix(tuple, i) = toFixed(vecValues[i]);
		}
		vecMatrix(tuple,2) = toFixed(1);
	}

	// //TODO: HACK 
	// if (pertubationID != -1)
	// {	
	// 	for (vtkIdType i = 0; i < vectors->GetNumberOfComponents() - 1; i++)
	// 		vecMatrix(pertubationID, i) = toFixed(0);
	// }
	
	//std::cout << " \t ######################################################################### " << std::endl;
	//std::cout << " \t\t" << vecMatrix(0, 0) << " " << vecMatrix(0, 1) << " " << vecMatrix(0, 2) << std::endl;
	//std::cout << " \t\t" << vecMatrix(1, 0) << " " << vecMatrix(1, 1) << " " << vecMatrix(1, 2) << std::endl;
	//std::cout << " \t\t" << vecMatrix(2, 0) << " " << vecMatrix(2, 1) << " " << vecMatrix(2, 2) << std::endl;
	//std::cout << " \t ######################################################################### " << std::endl;

	// 2. compute determinant sign
	double det = vecMatrix.determinant();

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		det *= -1;
		//if(det > 0)
		//	det *= -1;
	}
	//std::cout << "\t\t Determinat: " << det << " Swaps:" << swapOperations << std::endl;
	return det;
}

double CriticalPointExtractor::toFixed(double val)
{
	//TODO: Some magic here
	double ret = val * pow(10, 14);
	//auto tofixed = cnl::fixed_point<unsigned long, -10, 2>(val);
	//double ret = to_rep(tofixed);
	// static_assert(tofixed==val, "fixed-point type was unable to store the value");
	//std::cout << val << " " << to_rep(tofixed) << std::endl;

	return ret;
}

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

bool CriticalPointExtractor::DeterminatCounterClockWise(double det)
{
	if (det > 0)
		return true;
	else
		return false;
}
