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

  //Compute critical points
  double singularity[3] = {0.00, 0.00, 0.00};  

  auto start = std::chrono::steady_clock::now();
  CriticalPointExtractor cp_extractor(input, false, singularity);

  auto end = std::chrono::steady_clock::now();  
   std::cout << "[CriticalPointExtractor::Constructor] Elapsed time in milliseconds : "
 	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
 	  << " ms" << std::endl;

  start = std::chrono::steady_clock::now();
  cp_extractor.ComputeCriticalCells(output);
  end = std::chrono::steady_clock::now();
  std::cout << "[identify_critical_points] Unstable critical points: " << output->GetNumberOfCells() << std::endl;
  std::cout << "[identify_critical_points] Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;

  // Local cleanup done by every worker
  start = std::chrono::steady_clock::now();
  cp_extractor.CleanDuplicates(output);
  
  //Collect all critical cells per MPI rank
  vtkSmartPointer < vtkAggregateDataSetFilter > reducedData = vtkSmartPointer < vtkAggregateDataSetFilter >::New();
  reducedData->SetInputData(output);
  reducedData->Update();
  output->ShallowCopy(reducedData->GetPolyDataOutput());

  /// to-DO: post-processing DUPLICATE CLEANUP
  ///  Global cleanup. Only the master has results to polish
  cp_extractor.CleanDuplicates(output);
  end = std::chrono::steady_clock::now();
  std::cout << "[duplicate_cleanup] Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;
  std::cout << "[duplicate_cleanup] Stable critical points: " << output->GetNumberOfCells() << std::endl;

  /// TO-DO: comparison between preprocessing perturbation vs post-processing cleanup

  return 1;
}

CriticalPointExtractor::CriticalPointExtractor(vtkSmartPointer<vtkDataSet> input,
											   bool pertubate,
											   double *currentSingularity)
{
	//Configure openmp
	omp_set_num_threads(numThreads);

	//Store singularity
	singularity[0] = currentSingularity[0];
	singularity[1] = currentSingularity[1];
	singularity[2] = currentSingularity[2];

	vtkIdType numPoints = input->GetNumberOfPoints();
	ZERO_ID = numPoints + 1;

	//Allocate memory
	vecPointCoordinates.resize(numPoints);
	vecVectors.resize(numPoints);
	
	//Store vectors and point coordinates for internal usage
	vtkSmartPointer<vtkDataArray> vectors = input->GetPointData()->GetVectors();	

	#pragma omp parallel for
	for(vtkIdType i=0; i < numPoints; i++) 
	{
		double* position = new double[3];
		double* vector = new double[3];

		vectors->GetTuple(i,vector);
		input->GetPoint(i, position);

		if(pertubate)
		{
			Pertubate(vector, i);
			Pertubate(singularity, ZERO_ID);
		}
		
		vecVectors[i] = vector;
		vecPointCoordinates[i] = position;
	}

	//Check for dataset dimension and configure the index in the matrix which will be set to 1
	double bounds[6];
	input->GetBounds(bounds);
	double xDim = fabs(bounds[1] - bounds[0]);
	double yDim = fabs(bounds[3] - bounds[2]);
	double zDim = fabs(bounds[5] - bounds[4]);	

	iExchangeIndex = 3; 					//3D dataset 
	if (xDim == 0.0) iExchangeIndex = 0; 	//2D dataset with yz
	if (yDim == 0.0) iExchangeIndex = 1; 	//2D dataset with xz
	if (zDim == 0.0) iExchangeIndex = 2; 	//2D dataset with xy

	std::cout << "[CriticalPointExtractor::identify_critical_points] Checking " << input->GetNumberOfCells() << " cells for critical points " << std::endl;

	//Vector of critical cell ids
	std::vector<CriticalPoint> vecCriticalCellIDs;

	//Configure for parallel independent processing
    std::vector<vtkSmartPointer<vtkGenericCell>> vecCellPerThread;  //Cell for each thread
	std::vector<std::vector<vtkIdType>> idsForCell;					//Cell point ids for each thread

	//Get the cell type once (needed for correct allocation)
	vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
	vtkIdType numCells = input->GetNumberOfCells();
	input->GetCell(0, cell);
    vtkIdType cellType = cell->GetCellType();

	for(int x=0; x < numThreads;++x)
	{
		vecCellPerThread.push_back(vtkSmartPointer<vtkGenericCell>::New());
		if(VTK_PIXEL == cellType || VTK_QUAD == cellType || VTK_TRIANGLE == cellType)
			idsForCell.push_back(std::vector<vtkIdType>(3));
		else
			idsForCell.push_back(std::vector<vtkIdType>(4));
	}

	//Allocate size for cells which depends on input cell type
	if(VTK_PIXEL == cellType || VTK_QUAD == cellType) vecCellIds.resize(numCells * 2);
	else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType) vecCellIds.resize(numCells * 5);
	else if (VTK_TRIANGLE == cellType) vecCellIds.resize(numCells);
	else if (VTK_TETRA == cellType) vecCellIds.resize(numCells);

	#pragma omp parallel for
	for (vtkIdType i = 0; i < numCells; i++) {
		//Local variables per thread
		int threadIdx = omp_get_thread_num();//Thread ID

		//get the current cell
		input->GetCell(i, vecCellPerThread[threadIdx]);

		//Get the associated point ids for the cell 
		vtkSmartPointer<vtkIdList> ids = vecCellPerThread[threadIdx]->GetPointIds();
	
		if (VTK_PIXEL == cellType || VTK_QUAD == cellType)
		{
			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(1), ids->GetId(2)};
			vecCellIds[i * 2] = idsForCell[threadIdx];

			idsForCell[threadIdx] = { ids->GetId(1) , ids->GetId(3), ids->GetId(2)};
			vecCellIds[i * 2 + 1] = idsForCell[threadIdx];
		}
		else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType)
		{
			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(1), ids->GetId(3), ids->GetId(5)};
			vecCellIds[i * 5] = idsForCell[threadIdx];

			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(3), ids->GetId(5), ids->GetId(6)};
			vecCellIds[i * 5 + 1] = idsForCell[threadIdx];
		
			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(2), ids->GetId(3), ids->GetId(6)};
			vecCellIds[i * 5 + 2] = idsForCell[threadIdx];

			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(4), ids->GetId(5), ids->GetId(6)};
			vecCellIds[i * 5 + 3] = idsForCell[threadIdx];

			idsForCell[threadIdx] = { ids->GetId(3) , ids->GetId(5), ids->GetId(6), ids->GetId(7)};
			vecCellIds[i * 5 + 4] = idsForCell[threadIdx];
		}
		else if (VTK_TRIANGLE == cellType)
		{
			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(1), ids->GetId(2)};
			vecCellIds[i] = idsForCell[threadIdx];
		}
		else if (VTK_TETRA == cellType)
		{
			idsForCell[threadIdx] = { ids->GetId(0) , ids->GetId(1), ids->GetId(2), ids->GetId(3)};
			vecCellIds[i] = idsForCell[threadIdx];
		}
		else {
			std::cout << "[CriticalPointExtractor] Error: unknown cell type " << std::endl;
				continue;
		}	
		
	}
	std::cout << "[CriticalPointExtractor::identify_critical_points] Checking " << vecCellIds.size() << " simplices for critical points " << std::endl;

}

void CriticalPointExtractor::Pertubate(double *values, vtkIdType id) {
	// perturbation function f(e,i,j) = e
	// e ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	double eps = 1;
	int delta = 4;

	int i = id+1; //TODO: should be point id of the simplex
	int j = 0 + 1;
	double perturbation = /*(id*5+0)*eps;*/std::pow(eps,std::pow(2,i*delta-j));
	std::cout << "Pertubation " << std::pow(2,i*delta-j) << std::endl;
	//auto tofixed = cnl::fixed_point<long, -14>(values[0]+perturbation);
	//values[0] = to_rep(tofixed);	
	values[0] += perturbation;

	j = 0 + 2;
	perturbation = /*(id*5+1)*eps;*/std::pow(eps,std::pow(2,i*delta-j));
	//tofixed = cnl::fixed_point<long, -14>(values[1]+perturbation);
	//values[1] = to_rep(tofixed);
	values[1] += perturbation;

	j = 0 + 3;
	perturbation = /*(id*5+2)*eps;*/std::pow(eps,std::pow(2,i*delta-j));
	//tofixed = cnl::fixed_point<long, -14>(values[2]+perturbation);
	//values[2] = to_rep(tofixed);
	values[2] += perturbation;
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::ComputeCriticalCells(vtkSmartPointer<vtkDataSet> output) 
{
	vtkSmartPointer<vtkPolyData> outputData = vtkPolyData::SafeDownCast(output);

	vtkIdType cells_num = vecCellIds.size();
	
	//Vector of critical cell ids
	std::vector<CriticalPoint> vecCriticalCellIDs;

	int matrixSize = 3;
	if(iExchangeIndex == 3)
		matrixSize = 4;

    std::vector<DynamicMatrix> vecMatrices;
	for(int x=0; x < numThreads;++x)
		vecMatrices.push_back(DynamicMatrix(matrixSize,matrixSize));

	//std::cout << "[CriticalPointExtractor::identify_critical_points] Matrix size(" << matrixSize << "," << matrixSize << ")"<< std::endl;
	//std::cout << "[CriticalPointExtractor::identify_critical_points] Exchange index: " << iExchangeIndex << std::endl;
	//std::cout << "[CriticalPointExtractor::identify_critical_points] Identifing critical cells "<< std::endl;

	//Check for every cell if a critical point (passed singularity as argument) exists
#pragma omp parallel 
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();	//Thread ID

		//Remove synchronization with nowait
		#pragma omp for nowait
		for (vtkIdType i = 0; i < vecCellIds.size(); i++) {
			//If the cell contains a the singularity add them to the output and we can break
			CriticalPointType ret = PointInCell(vecCellIds[i], vecMatrices[threadIdx]);
			if (ret != REGULAR) {
				CriticalPoint tmp(i,ret);
#pragma omp critical
				vecCriticalCellIDs.push_back(tmp);
			}
		}
	}
	
	vtkSmartPointer<vtkPoints> pointArray = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIntArray> singularityType = vtkSmartPointer<vtkIntArray>::New();
	singularityType->SetName("singularityType");

	//Prepare output data sequentially. Insert every critical cell to output
	for (auto cellID : vecCriticalCellIDs)
	{
		std::vector<vtkIdType> vecVertexIds = vecCellIds[cellID.id];
		vtkSmartPointer<vtkIdList> newPointIDs = vtkSmartPointer<vtkIdList>::New();	

		for (int index = 0; index < vecVertexIds.size(); index++)
		{
			singularityType->InsertNextTuple1(cellID.type);
			double pCoords[3];

			pCoords[0] = vecPointCoordinates[vecCellIds[cellID.id][index]][0];
			pCoords[1] = vecPointCoordinates[vecCellIds[cellID.id][index]][1];
			pCoords[2] = vecPointCoordinates[vecCellIds[cellID.id][index]][2];
			vtkIdType pointID;
			pointID = pointArray->InsertNextPoint(pCoords[0], pCoords[1], pCoords[2]);
			newPointIDs->InsertNextId(pointID);
		}
		cellArray->InsertNextCell(newPointIDs);
	}

	//Add points and cells to polydata
	outputData->SetPoints(pointArray); 
	outputData->SetPolys(cellArray);
	outputData->GetPointData()->AddArray(singularityType);

	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(outputData);
	cleaner->Update();
	outputData->ShallowCopy(cleaner->GetOutput());	
}

void CriticalPointExtractor::CleanDuplicates(vtkSmartPointer<vtkPolyData> output) {

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

CriticalPointExtractor::CriticalPointType CriticalPointExtractor::PointInCell(std::vector<vtkIdType> &ids, DynamicMatrix &vecMatrix) {
	// 1. compute the initial sign of the determinant of the cell
	double targetDeterminant = ComputeDeterminant(ids, vecMatrix, false, 0);
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
		tmpDeterminat   = ComputeDeterminant(ids, vecMatrix, false, i);
		tmpDirection    = DeterminatCounterClockWise(tmpDeterminat);
		
		// 2.3. check if it changes --> if so return false
		if (targetDirection != tmpDirection)
		{
			return REGULAR; // regular cell
		}
	}

	double initialDeterminant = ComputeDeterminant(ids, vecMatrix, true);	
	bool initialDirection     = DeterminatCounterClockWise(initialDeterminant);	
	
	if (initialDirection != targetDirection)
	{
		return SADDLE; // we found a saddle
	}

	return SINGULARITY; // the cell is critical, since the sign never change
}

double CriticalPointExtractor::ComputeDeterminant(
	std::vector<vtkIdType> tmpIds,
	DynamicMatrix &vecMatrix,
	bool usePoints,
	long pertubationID
){
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
		vtkIdType pointID = tmpIds[tuple];
		if (pointID != ZERO_ID)
		{
			if(!usePoints)
			{
				vecValues[0] = vecVectors[pointID][0];
				vecValues[1] = vecVectors[pointID][1];
				vecValues[2] = vecVectors[pointID][2];
			}
			else
			{
				vecValues[0] = vecPointCoordinates[pointID][0];
				vecValues[1] = vecPointCoordinates[pointID][1];
				vecValues[2] = vecPointCoordinates[pointID][2];
			}
		}
		else
		{
			vecValues[0] = singularity[0];
			vecValues[1] = singularity[1];
			vecValues[2] = singularity[2];
		}

		for (vtkIdType i = 0; i < vecMatrix.cols(); i++) {
			vecMatrix(tuple, i) = vecValues[i];
		}
		vecMatrix(tuple, iExchangeIndex) = 1;
	}

	/*
	std::cout << " \t ######################################################################### " << std::endl;
	std::cout << " \t\tVertex IDs ";
	for (int x = 0; x < tmpIds.size(); x++)
		std::cout << tmpIds[x] << " ";
	std::cout << std::endl;
	std::cout << " \t\tMatrix: " << std::endl;
	std::cout << " \t\t ";
	for (int x = 0; x < vecMatrix.rows(); x++)
	{
		for (int y = 0; y < vecMatrix.cols(); y++)
			std::cout << vecMatrix(x, y) << " ";
		std::cout << std::endl;
		std::cout << " \t\t ";
	}
	std::cout << std::endl;
	std::cout << " \t ######################################################################### " << std::endl;
	*/

	// 2. compute determinant sign
	double det = vecMatrix.determinant();
	//long det = IntegerDeterminant(vecMatrix, vecMatrix.rows());

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		det *= -1;
	}
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
