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

#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkGenericCell.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkVoxel.h>
#include <vtkPixel.h>
#include <vtkAggregateDataSetFilter.h>

#include <sstream>
#include <chrono>
#include <cmath>


#include <omp.h>
#include <random>

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

  auto start = std::chrono::steady_clock::now();
  cp_extractor.identify_critical_points(input, output, singularities);
  auto end = std::chrono::steady_clock::now();

  //Collect all critical cells per MPI rank
  vtkSmartPointer < vtkAggregateDataSetFilter > reducedData = vtkSmartPointer < vtkAggregateDataSetFilter >::New();
  reducedData->SetInputData(output);
  reducedData->Update();
  output->ShallowCopy(reducedData->GetPolyDataOutput());
  std::cout << "Critical points found: " << output->GetNumberOfCells() << std::endl;
  std::cout << "Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;
  /// to-do MPI implementation
  ///  

  return 1;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::identify_critical_points(	
	vtkSmartPointer<vtkDataSet> input,
	vtkSmartPointer<vtkDataSet> output,
	std::vector<double*> singlarities
) {

	vtkSmartPointer<vtkPolyData> outputData = vtkPolyData::SafeDownCast(output);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  
	vtkIdType cells_num = input->GetNumberOfCells();
	///std::cout << "Checking " << cells_num << " cells for critical points " << std::endl;

	//Set zero id to end 
	ZERO_ID = input->GetNumberOfPoints();

	//Check for dataset dimension and configure zero vector exchange index
	double bounds[6];
	input->GetBounds(bounds);
	double xDim = fabs(bounds[1] - bounds[0]);
	double yDim = fabs(bounds[3] - bounds[2]);
	double zDim = fabs(bounds[5] - bounds[4]);	
	zeroDim = 3; //3D dataset 
	if (xDim == 0.0) zeroDim = 0; //2D dataset with yz
	if (yDim == 0.0) zeroDim = 1; //2D dataset with xz
	if (zDim == 0.0) zeroDim = 2; //2D dataset with xy

	// std::cout << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << std::endl;
	// std::cout << xDim << " " << yDim << " " << zDim << std::endl;
	// std::cout << input->GetCell(0)->GetCellDimension() <<std::endl;

	//create an eigen matrix
	int columns = 3;
	if (zeroDim == 3) columns = 4;

	//Prepare for parallel computation
	int numThreads = 12;
	omp_set_num_threads(numThreads);
	//Private cell for every thread to work on
	std::vector<vtkSmartPointer<vtkGenericCell>> vecCellPerThread; 
	std::vector<DynamicMatrix> vecMatrixPerThread; 

	//Vector of critical cell ids
	std::vector<vtkIdType> vecCriticalCellIDs;

	int rows = 3;
	if(zeroDim == 3)
		rows = 4;

	DynamicMatrix tmpMatrix = DynamicMatrix(rows, columns);

	for(int x=0; x < numThreads; x++) {
		vecCellPerThread.push_back(vtkSmartPointer<vtkGenericCell>::New());
		vecMatrixPerThread.push_back(tmpMatrix);
	}

		//Check for every cell if a critical point (passed singularity as argument) exists
#pragma omp parallel 
	{
		//Variables per thread
		int threadIdx = omp_get_thread_num();	//Thread ID
		double currentSingularity[3];			//The current singularity value to check for
		vtkIdType cellType;						//Current cell type

//Remove synchronization with nowait
#pragma omp for nowait
		for (vtkIdType i = 0; i < cells_num; i++) {
			//get the cell i and its point ids and store per thread private for processing
			input->GetCell(i, vecCellPerThread[threadIdx]);
			vtkSmartPointer<vtkIdList> ids = vecCellPerThread[threadIdx]->GetPointIds();
			
			std::vector<vtkSmartPointer<vtkIdList>> vecCells; //Cells to be processed per thread
			cellType = vecCellPerThread[threadIdx]->GetCellType();

			if (VTK_QUAD == cellType)
			{
				////CREATE 2 TRIANGLES				
				vtkSmartPointer<vtkIdList> tri_ids1 = vtkSmartPointer<vtkIdList>::New();
				tri_ids1->SetNumberOfIds(3);
				tri_ids1->SetId(0, ids->GetId(0));
				tri_ids1->SetId(1, ids->GetId(1));
				tri_ids1->SetId(2, ids->GetId(3));
				vecCells.push_back(tri_ids1);
				
				vtkSmartPointer<vtkIdList> tri_ids2 = vtkSmartPointer<vtkIdList>::New();
				tri_ids2->SetNumberOfIds(3);
				tri_ids2->SetId(0, ids->GetId(1));
				tri_ids2->SetId(1, ids->GetId(2));
				tri_ids2->SetId(2, ids->GetId(3));
				vecCells.push_back(tri_ids2);
			}
			else if (VTK_PIXEL == cellType)
			{
				// //CREATE 2 TRIANGLES
				/// DO NOT USE VTKTRIANGLE <---- MUCH SLOWER
				vtkSmartPointer<vtkIdList> tri_ids1 = vtkSmartPointer<vtkIdList>::New();	
				tri_ids1->SetNumberOfIds(3);
				tri_ids1->SetId(0, ids->GetId(0));
				tri_ids1->SetId(1, ids->GetId(1));
				tri_ids1->SetId(2, ids->GetId(2));
				vecCells.push_back(tri_ids1);

				vtkSmartPointer<vtkIdList> tri_ids2 = vtkSmartPointer<vtkIdList>::New();
				tri_ids2->SetNumberOfIds(3);
				tri_ids2->SetId(0, ids->GetId(1));
				tri_ids2->SetId(1, ids->GetId(3));
				tri_ids2->SetId(2, ids->GetId(2));
				vecCells.push_back(tri_ids2);
			}
			else if (VTK_VOXEL == cellType)
			{
				//Create 5 tetra
				vtkSmartPointer<vtkIdList> tet_ids1 = vtkSmartPointer<vtkIdList>::New();	
				tet_ids1->SetNumberOfIds(4);
				tet_ids1->SetId(0, ids->GetId(0));
				tet_ids1->SetId(1, ids->GetId(1));
				tet_ids1->SetId(2, ids->GetId(3));
				tet_ids1->SetId(3, ids->GetId(5));
				vecCells.push_back(tet_ids1);

				vtkSmartPointer<vtkIdList> tet_ids2 = vtkSmartPointer<vtkIdList>::New();	
				tet_ids2->SetNumberOfIds(4);
				tet_ids2->SetId(0, ids->GetId(0));
				tet_ids2->SetId(1, ids->GetId(3));
				tet_ids2->SetId(2, ids->GetId(5));
				tet_ids2->SetId(3, ids->GetId(6));
				vecCells.push_back(tet_ids2);

				vtkSmartPointer<vtkIdList> tet_ids3 = vtkSmartPointer<vtkIdList>::New();	
				tet_ids3->SetNumberOfIds(4);
				tet_ids3->SetId(0, ids->GetId(0));
				tet_ids3->SetId(1, ids->GetId(2));
				tet_ids3->SetId(2, ids->GetId(3));
				tet_ids3->SetId(3, ids->GetId(6));
				vecCells.push_back(tet_ids3);

				vtkSmartPointer<vtkIdList> tet_ids4 = vtkSmartPointer<vtkIdList>::New();	
				tet_ids4->SetNumberOfIds(4);
				tet_ids4->SetId(0, ids->GetId(0));
				tet_ids4->SetId(1, ids->GetId(4));
				tet_ids4->SetId(2, ids->GetId(5));
				tet_ids4->SetId(3, ids->GetId(6));
				vecCells.push_back(tet_ids4);

				vtkSmartPointer<vtkIdList> tet_ids5 = vtkSmartPointer<vtkIdList>::New();	
				tet_ids5->SetNumberOfIds(4);
				tet_ids5->SetId(0, ids->GetId(3));
				tet_ids5->SetId(1, ids->GetId(5));
				tet_ids5->SetId(2, ids->GetId(6));
				tet_ids5->SetId(3, ids->GetId(7));
				vecCells.push_back(tet_ids5);
			}
			else if (VTK_TRIANGLE == cellType)
			{
				vtkSmartPointer<vtkIdList> tri_ids1 = vtkSmartPointer<vtkIdList>::New();	
				tri_ids1->SetNumberOfIds(3);
				tri_ids1->SetId(0, ids->GetId(0));
				tri_ids1->SetId(1, ids->GetId(1));
				tri_ids1->SetId(2, ids->GetId(2));
				vecCells.push_back(tri_ids1);
			}
			else if (VTK_TETRA == cellType)
			{
				vtkSmartPointer<vtkIdList> tet = vtkSmartPointer<vtkIdList>::New();	
				tet->SetNumberOfIds(4);
				tet->SetId(0, ids->GetId(0));
				tet->SetId(1, ids->GetId(1));
				tet->SetId(2, ids->GetId(2));
				tet->SetId(3, ids->GetId(3));
				vecCells.push_back(tet);
			}
			else {
				std::cout << "[CriticalPointExtractor] Error: unknown cell type. Number of point ids " << ids->GetNumberOfIds() << std::endl;
			}

			//Compute if one of the cells contains the given singularity (normally 0 in any dimension)
			for (auto cellFromVec : vecCells)
			{
				currentSingularity[0] = 0;
				currentSingularity[1] = 0;
				currentSingularity[2] = 0;

				//If the cell contains a the singularity add them to the output and we can break
				if (PointInCell(cellFromVec, input, currentSingularity, vecMatrixPerThread[threadIdx])) {
#pragma omp critical
					vecCriticalCellIDs.push_back(i);
					break;
				}
			}
		}
	}
	//Prepare output data sequentially. Insert every critical cell to output
	for (auto cellID : vecCriticalCellIDs)
	{
		vtkSmartPointer<vtkCell> cell = input->GetCell(cellID);
		vtkSmartPointer<vtkIdList> ids = cell->GetPointIds();
		vtkSmartPointer<vtkIdList> new_ids = vtkSmartPointer<vtkIdList>::New();;
		for (int index = 0; index < ids->GetNumberOfIds(); index++)
		{
			double pCoords[3];
			vtkIdType newCellID;
			input->GetPoint(ids->GetId(index), pCoords);
			newCellID = points->InsertNextPoint(pCoords[0], pCoords[1], pCoords[2]);
			new_ids->InsertNextId(newCellID);
		}
		cells->InsertNextCell(new_ids);
	}
	//Add points and cells to polydata
	outputData->SetPoints(points); 
	outputData->SetPolys(cells);

	//std::cout << "Critical points found: " << outputData->GetNumberOfCells() << std::endl;
}


bool CriticalPointExtractor::PointInCell(vtkSmartPointer<vtkIdList> ids, vtkSmartPointer<vtkDataSet> grid, double* currentSingularity, DynamicMatrix &vecMatrix) {
	//std::cout << " ################### Point in Cell ###################################################### " << std::endl;
	// 1. compute the initial sign of the determinant of the cell
	double initialDeterminant = Positive(ids, grid, currentSingularity, vecMatrix);	
	bool initialDirection     = DeterminatCounterClockWise(initialDeterminant);	

	//Check for non data values (vector is zero and determinant also) 
	if (initialDeterminant == 0)
	{
		return false;
	}

	double tmpDeterminat;
	bool tmpDirection;
	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0) given as i to Positive function
	for (int i = 0; i < ids->GetNumberOfIds(); i++) {
		// 2.2. compute the determinant sign again 
		tmpDeterminat   = Positive(ids, grid, currentSingularity, vecMatrix, i);
		tmpDirection    = DeterminatCounterClockWise(tmpDeterminat);
		
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
double CriticalPointExtractor::Positive(
	vtkSmartPointer<vtkIdList> ids,
	vtkSmartPointer<vtkDataSet> grid,
	double currentSingularity[3],
	DynamicMatrix &vecMatrix,
	long pertubationID
){
	//Copy ids for local modification
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();
	std::vector<vtkIdType> tmpIds(ids->GetNumberOfIds());
	for (vtkIdType tuple = 0; tuple < ids->GetNumberOfIds(); tuple++)
		tmpIds[tuple] = ids->GetId(tuple);

	//Exchanges every facet with the zero vector
	if (pertubationID != -1)
	{
		tmpIds[pertubationID] = ZERO_ID;
	}

	// 1. Sort and check swap operations (check)
    // TODO: More generic version required. How to handle per pertubation
	int swapOperations =  Sort(tmpIds);
	
	double vecValues[3];
	for (std::size_t tuple = 0; tuple < tmpIds.size(); tuple++) 
	{
		if (tmpIds[tuple] != ZERO_ID)
		{
			vectors->GetTuple(tmpIds[tuple], vecValues);
		}
		else
		{
			vecValues[0] = currentSingularity[0];
			vecValues[1] = currentSingularity[1];
			vecValues[2] = currentSingularity[2];
		}

		for (vtkIdType i = 0; i < vecMatrix.cols(); i++) {
			//TODO: long long to fixed precision
			vecMatrix(tuple, i) = toFixed(vecValues[i]);
		}
		vecMatrix(tuple, zeroDim) = toFixed(tmp);
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
	std::cout << " \t ######################################################################### " << std::endl;*/

	// 2. compute determinant sign
	double det = vecMatrix.determinant();

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		det *= -1;
	}
	//std::cout << "\t\t Determinat: " << det << " Swaps:" << swapOperations << std::endl;
	return det;
}

double CriticalPointExtractor::toFixed(double& val)
{
	//TODO: Some magic here
	//long long ret =  val * pow(10, 14);
	double ret = val;

	return ret;
}

int  CriticalPointExtractor::Sort(std::vector<vtkIdType> &ids)
{
	if (ids.size() == 3) //Triangle
	{
		return Sort3(ids);
	}
	else if (ids.size() == 4) //TETRA 
	{
		return Sort4(ids);
	}
	else
	{
		std::cout << "Warning cell type currently not supported" << std::endl;
		return 0;
	}
}

int  CriticalPointExtractor::Sort3(std::vector<vtkIdType> &ids)
{
	unsigned int swaps = 0;
	if (ids[0] > ids[1])
	{
		std::swap(ids[0],ids[1]);
		swaps++;
	}

	if (ids[1] > ids[2])
	{
		std::swap(ids[1],ids[2]);
		swaps++;

		if (ids[0] > ids[1])
		{
			std::swap(ids[0],ids[1]);
			swaps++;
		}
	}
	return swaps;
}

int  CriticalPointExtractor::Sort4(std::vector<vtkIdType> &ids)
{
	unsigned int swaps = 0;

	if (ids[0] > ids[1])
	{
		std::swap(ids[0],ids[1]);
		swaps++;
	}

	if (ids[1] > ids[2])
	{
		std::swap(ids[1],ids[2]);
		swaps++;

		if (ids[0] > ids[1])
		{
			std::swap(ids[0],ids[1]);
			swaps++;
		}
	}

	if (ids[3] < ids[2])
	{
		if (ids[3] < ids[0])
		{
			std::swap(ids[2],ids[3]);
			swaps++;

			std::swap(ids[1],ids[2]);
			swaps++;

			std::swap(ids[0],ids[1]);
			swaps++;
		}
		else if (ids[3] < ids[1])
		{
			std::swap(ids[2],ids[3]);
			swaps++;

			std::swap(ids[1],ids[2]);
			swaps++;
		}
		else
		{
			std::swap(ids[2],ids[3]);
			swaps++;
		}
	}
	return swaps;
}

bool CriticalPointExtractor::DeterminatCounterClockWise(double& det)
{
	if (det > 0)
		return true;
	else
		return false;
}
