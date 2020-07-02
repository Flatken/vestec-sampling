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
#include <vtkCleanPolyData.h>

#include <sstream>
#include <chrono>
#include <cmath>


#include <omp.h>
#include <random>
#include <array>

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

  CriticalPointExtractor cp_extractor;

  //Set zero id to end 
  cp_extractor.ZERO_ID = input->GetNumberOfPoints() + 1;

  //Compute critical points
  double singularity[3] = {0.00, 0.00, 0.00};  

  auto start = std::chrono::steady_clock::now();
//   cp_extractor.toFixed(singularity,cp_extractor.ZERO_ID);

  //cp_extractor.toFixed(cp_extractor.ONE,cp_extractor.ZERO_ID+1);

  /// TO-DO: GLOBAL PERTURBATION ON THE DATA
  /// ----> DOUBLE TO FIXED PRECISION CONVERSION 
  /// THIS CAN BE DISABLE
//   cp_extractor.perturbate(input);
  auto end = std::chrono::steady_clock::now();  
//   std::cout << "[perturbate] Elapsed time in milliseconds : "
// 	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
// 	  << " ms" << std::endl;

  start = std::chrono::steady_clock::now();
  cp_extractor.identify_critical_points(input, output, singularity);
  end = std::chrono::steady_clock::now();
  std::cout << "[identify_critical_points] Unstable critical points: " << output->GetNumberOfCells() << std::endl;
  std::cout << "[identify_critical_points] Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;

  // Local cleanup done by every worker
  start = std::chrono::steady_clock::now();
  cp_extractor.duplicate_cleanup(output);
  
  //Collect all critical cells per MPI rank
  vtkSmartPointer < vtkAggregateDataSetFilter > reducedData = vtkSmartPointer < vtkAggregateDataSetFilter >::New();
  reducedData->SetInputData(output);
  reducedData->Update();
  output->ShallowCopy(reducedData->GetPolyDataOutput());

  /// to-DO: post-processing DUPLICATE CLEANUP
  ///  Global cleanup. Only the master has results to polish
  cp_extractor.duplicate_cleanup(output);
  end = std::chrono::steady_clock::now();
  std::cout << "[duplicate_cleanup] Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;
  std::cout << "[duplicate_cleanup] Stable critical points: " << output->GetNumberOfCells() << std::endl;

  /// TO-DO: comparison between preprocessing perturbation vs post-processing cleanup

  return 1;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::perturbate(vtkSmartPointer<vtkDataSet> grid) {
		//Copy ids for local modification
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();	
	for(vtkIdType i=0; i<vectors->GetNumberOfTuples(); i++) {
		double values[3];		
		vectors->GetTuple(i,values);	
		toFixed(values,i);
		vectors->SetTuple(i,values);
	}
}

void CriticalPointExtractor::toFixed(double *values, vtkIdType id) {
	// longi := int (value * 10^a).
	/*char str[50];
		
	long long w = 15;
	long long a = 14;

	sprintf(str,"%0.*f",a,values[0]);
	values[0] = atof(str) * (double) std::pow(10,a);
	sprintf(str, "%0.1f", floor (values[0]));
	values[0] = (long long) atof(str);
	
	sprintf(str,"%0.*f",a,values[1]);
	values[1] = atof(str) * (double) std::pow(10,a);
	sprintf(str, "%0.1f", floor (values[1]));
	values[1] = (long long) atof(str);
	
	sprintf(str,"%0.*f",a,values[2]);
	values[2] = atof(str) * (double) std::pow(10,a);
	sprintf(str, "%0.1f", floor (values[2]));
	values[2] = (long long) atof(str);*/	

	// perturbation function f(e,i,j) = e
	// e ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	double eps = 0.05;
	int delta = 5;

	int i = id+1;
	int j = 0 + 1;
	double perturbation = /*(id*5+0)*eps;*/std::pow(eps,std::pow(2,i*delta-j));
	// values[0] *= std::pow(10,14);
	// auto tofixed = cnl::fixed_point<long long, -14>(values[0]+perturbation);
	// values[0] = to_rep(tofixed);	
	values[0] += perturbation;

	j = 0 + 2;
	perturbation = /*(id*5+1)*eps;*/std::pow(eps,std::pow(2,i*delta-j));
	// values[1] *= std::pow(10,14);
	// tofixed = cnl::fixed_point<long long, -14>(values[1]+perturbation);
	// values[1] = to_rep(tofixed);
	values[1] += perturbation;

	j = 0 + 3;
	perturbation = /*(id*5+2)*eps;*/std::pow(eps,std::pow(2,i*delta-j));
	// values[2] *= std::pow(10,14);
	// tofixed = cnl::fixed_point<long long, -14>(values[2]+perturbation);
	// values[2] = to_rep(tofixed);
	values[2] += perturbation;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::identify_critical_points(	
	vtkSmartPointer<vtkDataSet> input,
	vtkSmartPointer<vtkDataSet> output,
	double* singlarity
) {

	vtkSmartPointer<vtkPolyData> outputData = vtkPolyData::SafeDownCast(output);
	  
	vtkIdType cells_num = input->GetNumberOfCells();
	std::cout << "[CriticalPointExtractor::identify_critical_points] Checking " << cells_num << " cells for critical points " << std::endl;

	

	//Check for dataset dimension and configure zero vector exchange index
	double bounds[6];
	input->GetBounds(bounds);
	double xDim = fabs(bounds[1] - bounds[0]);
	double yDim = fabs(bounds[3] - bounds[2]);
	double zDim = fabs(bounds[5] - bounds[4]);	

	iExchangeIndex = 3; //3D dataset 

	if (xDim == 0.0) iExchangeIndex = 0; //2D dataset with yz
	if (yDim == 0.0) iExchangeIndex = 1; //2D dataset with xz
	if (zDim == 0.0) iExchangeIndex = 2; //2D dataset with xy

	//create an eigen matrix
	int iMatrixColumns = 3;
	if (iExchangeIndex == 3) iMatrixColumns = 4;

	//Prepare for parallel computation
	int numThreads = 12;
	omp_set_num_threads(numThreads);

	//Create a thread private cell for concurrent computation
	std::vector<vtkGenericCell*> vecCellPerThread; 
	std::vector<DynamicMatrix> vecMatrixPerThread; 

	//Vector of critical cell ids
	std::vector<CriticalPoint> vecCriticalCellIDs;

	int rows = 3;
	if(iExchangeIndex == 3)
		rows = 4;

	std::cout << "[CriticalPointExtractor::identify_critical_points] Matrix size(" << rows << "," << iMatrixColumns << ")"<< std::endl;
	std::cout << "[CriticalPointExtractor::identify_critical_points] Exchange index: " << iExchangeIndex << std::endl;
	for(int x=0; x < numThreads; x++) {
		vecCellPerThread.push_back(vtkGenericCell::New());
		vecMatrixPerThread.push_back(DynamicMatrix(rows, iMatrixColumns));
	}

	//Need to initialize for parallel processing: so call from main thread one time
	input->GetCell(0, vecCellPerThread[0]);

	std::cout << "[CriticalPointExtractor::identify_critical_points] Identifing critical cells "<< std::endl;
	//Check for every cell if a critical point (passed singularity as argument) exists
#pragma omp parallel 
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();	//Thread ID
		double currentSingularity[3];			//The current singularity value to check for
		vtkIdType cellType = -1;				//Current cell type
		std::array<std::vector<vtkIdType>, 5> arrayCells{ std::vector<vtkIdType>(3), std::vector<vtkIdType>(3),
														  std::vector<vtkIdType>(3), std::vector<vtkIdType>(3),
														  std::vector<vtkIdType>(3) };
		int generatedCells = 1;

		currentSingularity[0] = singlarity[0];
		currentSingularity[1] = singlarity[1];
		currentSingularity[2] = singlarity[2];

		//Remove synchronization with nowait
		#pragma omp for nowait
		for (vtkIdType i = 0; i < cells_num; i++) {
			//get the current cell
			input->GetCell(i, vecCellPerThread[threadIdx]);
			
			//Get the associated point ids for the cell 
			vtkSmartPointer<vtkIdList> ids = vecCellPerThread[threadIdx]->GetPointIds();
		
			//Get the cell type
			cellType = vecCellPerThread[threadIdx]->GetCellType();
			
			if (VTK_PIXEL == cellType || VTK_QUAD == cellType)
			{
				generatedCells = 2;
				arrayCells[0] = {ids->GetId(0) , ids->GetId(1), ids->GetId(2)};
				arrayCells[1] = {ids->GetId(1) , ids->GetId(3), ids->GetId(2)};
			}
			else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType)
			{
				//Create 5 tetra
				generatedCells = 5;

				arrayCells[0].resize(4);
				arrayCells[0][0] = ids->GetId(0);
				arrayCells[0][1] = ids->GetId(1);
				arrayCells[0][2] = ids->GetId(3);
				arrayCells[0][3] = ids->GetId(5);
				
				arrayCells[1].resize(4);
				arrayCells[1][0] = ids->GetId(0);
				arrayCells[1][1] = ids->GetId(3);
				arrayCells[1][2] = ids->GetId(5);
				arrayCells[1][3] = ids->GetId(6);
				
				arrayCells[2].resize(4);
				arrayCells[2][0] = ids->GetId(0);
				arrayCells[2][1] = ids->GetId(2);
				arrayCells[2][2] = ids->GetId(3);
				arrayCells[2][3] = ids->GetId(6);
				
				arrayCells[3].resize(4);
				arrayCells[3][0] = ids->GetId(0);
				arrayCells[3][1] = ids->GetId(4);
				arrayCells[3][2] = ids->GetId(5);
				arrayCells[3][3] = ids->GetId(6);
				
				arrayCells[4].resize(4);
				arrayCells[4][0] = ids->GetId(3);
				arrayCells[4][1] = ids->GetId(5);
				arrayCells[4][2] = ids->GetId(6);
				arrayCells[4][3] = ids->GetId(7);
			}
			else if (VTK_TRIANGLE == cellType)
			{
				generatedCells = 1;
				arrayCells[0] = {ids->GetId(0), ids->GetId(1), ids->GetId(2)};
			}
			else if (VTK_TETRA == cellType)
			{
				generatedCells = 1;
				arrayCells[0] = {ids->GetId(0), ids->GetId(1), ids->GetId(2), ids->GetId(3)};
			}
			else {
				std::cout << "[CriticalPointExtractor] Error: unknown cell type " << std::endl;
				continue;
			}

			//Compute if one of the cells contains the given singularity (normally 0 in any dimension)
			for(vtkIdType nb = 0; nb < generatedCells; nb++)
			{
				//If the cell contains a the singularity add them to the output and we can break
				CriticalPointType ret = PointInCell(arrayCells[nb], input, currentSingularity, vecMatrixPerThread[threadIdx]);
				if (ret != REGULAR) {
					CriticalPoint tmp(i,ret);
#pragma omp critical
					vecCriticalCellIDs.push_back(tmp);
					break;
				}
			}
		}
	}
	std::cout << "[CriticalPointExtractor::identify_critical_points] Identifing critical cells done "<< std::endl;
	std::cout << "[CriticalPointExtractor::identify_critical_points] Now preparing output data " << std::endl;

	vtkSmartPointer<vtkPoints> pointArray = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIntArray> singularityType = vtkSmartPointer<vtkIntArray>::New();
	singularityType->SetName("singularityType");

	//Prepare output data sequentially. Insert every critical cell to output
	for (auto cellID : vecCriticalCellIDs)
	{
		vtkSmartPointer<vtkCell> cell  = input->GetCell(cellID.id);
		vtkSmartPointer<vtkIdList> oldPointIDs = cell->GetPointIds();
		vtkSmartPointer<vtkIdList> newPointIDs = vtkSmartPointer<vtkIdList>::New();	

		for (int index = 0; index < oldPointIDs->GetNumberOfIds(); index++)
		{
			singularityType->InsertNextTuple1(cellID.type);
			double pCoords[3];
			vtkIdType pointID;
			input->GetPoint(oldPointIDs->GetId(index), pCoords);
			pointID = pointArray->InsertNextPoint(pCoords[0], pCoords[1], pCoords[2]);
			newPointIDs->InsertNextId(pointID);
		}
		cellArray->InsertNextCell(newPointIDs);
		//outputData->InsertNextCell(newPointIDs);
	}
	//Add points and cells to polydata
	outputData->SetPoints(pointArray); 
	outputData->SetPolys(cellArray);
	outputData->GetPointData()->AddArray(singularityType);

	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(outputData);
	cleaner->Update();
	outputData->ShallowCopy(cleaner->GetOutput());	
	//std::cout << "Critical points found: " << outputData->GetNumberOfCells() << std::endl;
}

void CriticalPointExtractor::duplicate_cleanup(vtkSmartPointer<vtkPolyData> output) {

	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < output->GetNumberOfPoints(); i++) {
		output->GetPointCells(i,ids);
		if(ids->GetNumberOfIds() > 1){
			for (vtkIdType id=1; id < ids->GetNumberOfIds(); id++) {			
				output->DeleteCell(ids->GetId(id));
			}
		}		
}
	output->RemoveDeletedCells();	

	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->PointMergingOff();
	cleaner->SetInputData(output);
	cleaner->Update();
	output->ShallowCopy(cleaner->GetOutput());		
}

CriticalPointExtractor::CriticalPointType CriticalPointExtractor::PointInCell(std::vector<vtkIdType> &ids, vtkSmartPointer<vtkDataSet> grid, double* currentSingularity, DynamicMatrix &vecMatrix) {
	//std::cout << " ################### Point in Cell ###################################################### " << std::endl;
	// 1. compute the initial sign of the determinant of the cell
	double targetDeterminant = ComputeDeterminant(ids, grid, currentSingularity, vecMatrix, false, 0);
	bool targetDirection = DeterminatCounterClockWise(targetDeterminant);
	//Check for non data values (vector is zero and determinant also) 
	if (targetDeterminant == 0)
	{
		return REGULAR;
	}

	double tmpDeterminat;
	bool tmpDirection;
	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0) given as i to Positive function
	for (int i = 1; i < ids.size(); i++) {
		// 2.2. compute the determinant sign again 
		tmpDeterminat   = ComputeDeterminant(ids, grid, currentSingularity, vecMatrix, false, i);
		tmpDirection    = DeterminatCounterClockWise(tmpDeterminat);
		
		// 2.3. check if it changes --> if so return false
		if (targetDirection != tmpDirection)
		{
			return REGULAR; // regular cell
		}
	}

	double initialDeterminant = ComputeDeterminant(ids, grid, currentSingularity, vecMatrix, true);	
	bool initialDirection     = DeterminatCounterClockWise(initialDeterminant);	
	// //Check for non data values (vector is zero and determinant also) 
	// if (initialDeterminant == 0)
	// {
	// 	return REGULAR;
	// }
	
	//
	if (initialDirection != targetDirection)
	{
		return SADDLE; // we found a saddle
	}

	//std::cout << " \t\t Cell is critical "  << std::endl;
	return SINGULARITY; // the cell is critical, since the sign never change
}

double CriticalPointExtractor::ComputeDeterminant(
	std::vector<vtkIdType> tmpIds,
	vtkSmartPointer<vtkDataSet> grid,
	double currentSingularity[3],
	DynamicMatrix &vecMatrix,
	bool usePoints,
	long pertubationID
){
	//Copy ids for local modification
	vtkSmartPointer<vtkDataArray> vectors = grid->GetPointData()->GetVectors();		

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
			if(!usePoints)
				vectors->GetTuple(tmpIds[tuple], vecValues);
			else
				grid->GetPoint(tmpIds[tuple], vecValues);
		}
		else
		{
			vecValues[0] = currentSingularity[0];
			vecValues[1] = currentSingularity[1];
			vecValues[2] = currentSingularity[2];
		}

		for (vtkIdType i = 0; i < vecMatrix.cols(); i++) {
			//TODO: long long to fixed precision
			vecMatrix(tuple, i) = vecValues[i];
		}
		vecMatrix(tuple, iExchangeIndex) = ONE[0];
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

int CriticalPointExtractor::Sort(std::vector<vtkIdType> &ids)
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

int CriticalPointExtractor::Sort3(std::vector<vtkIdType> &ids)
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

int CriticalPointExtractor::Sort4(std::vector<vtkIdType> &ids)
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
