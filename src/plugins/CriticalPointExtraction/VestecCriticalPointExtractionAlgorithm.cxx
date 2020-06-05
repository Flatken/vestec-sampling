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
  double p1[3] = {0.00, 0.00, 0};
  double p2[3] = {0.0, 0.001, 0};
  double p3[3] = {0.001, 0.0, 0};

  singularities.push_back(p1);
  //singularities.push_back(p2);
  //singularities.push_back(p3);

  CriticalPointExtractor cp_extractor;
  cp_extractor.identify_critical_points(input, output, singularities);
  /// to-do MPI implementation
  ///  

  return 1;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::identify_critical_points(	vtkSmartPointer<vtkDataSet> input,
																				vtkSmartPointer<vtkDataSet> output, std::vector<double*> singlarities) {

	int cp = 0;
	vtkSmartPointer<vtkPolyData> outputData = vtkPolyData::SafeDownCast(output);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  
	vtkIdType cells_num = input->GetNumberOfCells();
	std::cout << "Checking " << cells_num << " cells for critical points " << std::endl;

	//Set zero id to end 
	ZERO_ID = input->GetNumberOfPoints();

	//Prepare for parallel computation
	int numThreads = 12;
	omp_set_num_threads(numThreads);
	std::vector<vtkSmartPointer<vtkGenericCell>> vecCellPerThread;
	for(int x=0; x < numThreads; x++)
		vecCellPerThread.push_back(vtkSmartPointer<vtkGenericCell>::New());
#

		//Check for every cell if a critical point (passed singularity as argument) exists
//#pragma omp parallel for
		for (vtkIdType i = 0; i < cells_num; i++) {
			//get the cell for the current thread
			int threadIdx = omp_get_thread_num();
			input->GetCell(i, vecCellPerThread[threadIdx]);

			//Compute if cell contains the given singularity (normally 0 in any dimension)
			for (auto singularity : singlarities)
			{
				double currentSingularity[3];
				currentSingularity[0] = singularity[0];
				currentSingularity[1] = singularity[1];
				currentSingularity[2] = singularity[2];

				//If the cell contains a the singularity add them to the output
				if (PointInCell(vecCellPerThread[threadIdx], input, currentSingularity)) {
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
					cp++;
				}
				//char n;
				//std::cout << "Press for next";
				//std::cin >> n;

			}
		}
	

  //Add points and cells to polydata
  outputData->SetPoints(points); 
  outputData->SetPolys(cells);
  std::cout << "Critical points found: " << cp << std::endl;
}

bool CriticalPointExtractor::PointInCell(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid, double* currentSingularity) {
	vtkSmartPointer<vtkIdList> ids = cell->GetPointIds();

	//std::cout << " ################### Point in Cell ###################################################### " << std::endl;
	// 1. compute the initial sign of the determinant of the cell
	double initialDeterminant = Positive(ids, grid, currentSingularity);
	bool initialDirection     = DeterminatCounterClockWise(initialDeterminant);

	//Check for non data values (vector is zero and determinant also) 
	if (initialDeterminant == 0)
	{
		return false;
	}

	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0) given as i to Positive function
	for (int i = 0; i < ids->GetNumberOfIds(); i++) {
		// 2.2. compute the determinant sign again 
		double tmpDeterminat = Positive(ids, grid, currentSingularity, i);
		bool tmpDirection    = DeterminatCounterClockWise(tmpDeterminat);

		// 2.3. check if it changes --> if so return false
		if (initialDirection != tmpDirection)
		{
			return false;
		}
	}
	//std::cout << " \t\t Cell is critical "  << std::endl;
	return true; // the cell is critical, since the sign never change
}

/// can we pass to Positive directly the determinant matrix instead of the cell?
double CriticalPointExtractor::Positive(vtkSmartPointer<vtkIdList> ids, vtkSmartPointer<vtkDataSet> grid, double currentSingularity[3], long pertubationID){
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();
	vtkSmartPointer<vtkIdList> tmpIds = vtkSmartPointer<vtkIdList>::New();
	tmpIds->DeepCopy(ids);

	//Check for dataset dimension and configure zero vector exchange index
	double bounds[6];
	grid->GetBounds(bounds);
	double xDim = fabs(bounds[1] - bounds[1]);
	double yDim = fabs(bounds[3] - bounds[2]);
	double zDim = fabs(bounds[5] - bounds[4]);
	int zeroDim = 3; //3D dataset 
	if (xDim == 0) zeroDim = 0; //2D dataset with yz
	if (yDim == 0) zeroDim = 1; //2D dataset with xz
	if (zDim == 0) zeroDim = 2; //2D dataset with xy

	//Exchanges every facet with the zero vector
	if (pertubationID != -1)
	{
		tmpIds->SetId(pertubationID, ZERO_ID);
	}

	// 1. Sort and check swap operations (check)
    // TODO: More generic version required. How to handle per pertubation
	int swapOperations = Sort(tmpIds);

	//create an eigen matrix
	int columns = 3;
	if (zeroDim == 3) columns = 4;

	DynamicMatrix vecMatrix(tmpIds->GetNumberOfIds(), columns);
	
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
			vecMatrix(tuple, i) = toFixed(vecValues[i]);
		}
		vecMatrix(tuple, zeroDim) = toFixed(1);
	}

	/*std::cout << " \t ######################################################################### " << std::endl;
	std::cout << " \t\tVertex IDs ";
	for (int x = 0; x < tmpIds->GetNumberOfIds(); x++)
		std::cout << tmpIds->GetId(x) << " ";
	std::cout << std::endl;
	std::cout << " \t\tMatrix: " << std::endl;
	std::cout << " \t\t ";
	for (int x = 0; x < tmpIds->GetNumberOfIds(); x++)
	{
		for (int y = 0; y < columns; y++)
			std::cout << vecMatrix(x, y) << " ";
		std::cout << std::endl;
		std::cout << " \t\t ";
	}
	std::cout << std::endl;
	std::cout << " \t ######################################################################### " << std::endl;
*/
	// 2. compute determinant sign
	double det = vecMatrix.determinant();

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		det *= -1;

		//if(det > 0)
			//det *= -1;
	}
	//std::cout << "\t\t Determinat: " << det << " Swaps:" << swapOperations << std::endl;
	return det;
}

double CriticalPointExtractor::toFixed(double val)
{
	//TODO: Some magic here
	//long long ret =  val * pow(10, 14);
	double ret = val;

	return ret;
}

int CriticalPointExtractor::Sort(vtkSmartPointer<vtkIdList> ids)
{
	if (ids->GetNumberOfIds() == 3) //Triangle
	{
		return Sort3(ids);
	}
	else if (ids->GetNumberOfIds() == 4) //QUAD or TERAHERDON
	{
		return Sort4(ids);
	}
	else
	{
		std::cout << "Warning cell type currently not supported" << std::endl;
		return 0;
	}
}

int CriticalPointExtractor::Sort3(vtkSmartPointer<vtkIdList> ids)
{
	vtkIdType tmp;
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

int  CriticalPointExtractor::Sort4(vtkSmartPointer<vtkIdList> ids)
{
	unsigned int swaps = 0;
	vtkIdType tmp;

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

	if (ids->GetId(3) < ids->GetId(2))
	{
		if (ids->GetId(3) < ids->GetId(0))
		{
			tmp = ids->GetId(2);
			ids->SetId(2, ids->GetId(3));
			ids->SetId(3, tmp);
			swaps++;

			tmp = ids->GetId(1);
			ids->SetId(1, ids->GetId(2));
			ids->SetId(2, tmp);
			swaps++;

			tmp = ids->GetId(0);
			ids->SetId(0, ids->GetId(1));
			ids->SetId(1, tmp);
			swaps++;
		}
		else if (ids->GetId(3) < ids->GetId(1))
		{
			tmp = ids->GetId(2);
			ids->SetId(2, ids->GetId(3));
			ids->SetId(3, tmp);
			swaps++;

			tmp = ids->GetId(1);
			ids->SetId(1, ids->GetId(2));
			ids->SetId(2, tmp);
			swaps++;
		}
		else
		{
			tmp = ids->GetId(2);
			ids->SetId(2, ids->GetId(3));
			ids->SetId(3, tmp);
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
