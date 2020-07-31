#include "VestecCriticalPointExtractionAlgorithm.h"

#include <vtkPoints.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCommand.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMultiProcessController.h>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkGenericCell.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>
#include <vtkVoxel.h>
#include <vtkPixel.h>
#include <vtkAggregateDataSetFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkCleanUnstructuredGrid.h>
#include <vtkDataSetWriter.h>

#include <sstream>
#include <chrono>
#include <cmath>
#include <thread>

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
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
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
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::GetData(outputVector, 0);

  //Set active array
  vtkDataArray *inScalars = this->GetInputArrayToProcess(0, inputVector);
  input->GetPointData()->SetActiveVectors(inScalars->GetName());

  //Get number of points from each MPI process
  //calculate local start index
  //calculate total number of vertices

  //Compute critical points
  double singularity[3] = {0.000, 0.000, 0.000};  

  auto start = std::chrono::steady_clock::now();
  CriticalPointExtractor cp_extractor(input, singularity); //perturbation is ON by default

  auto end = std::chrono::steady_clock::now();  
   std::cout << "[CriticalPointExtractor::Constructor] Elapsed time in milliseconds : "
 	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
 	  << " ms" << std::endl;

  start = std::chrono::steady_clock::now();
  cp_extractor.ComputeCriticalCells(output);
  end = std::chrono::steady_clock::now();
  std::cout << "[identify_critical_points] Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;
  std::cout << "[identify_critical_points] Unstable critical cells: " << output->GetNumberOfCells() << std::endl;

  // Local cleanup done by every worker
  start = std::chrono::steady_clock::now();
  
  //Collect all critical cells per MPI rank
  vtkSmartPointer < vtkAggregateDataSetFilter > reducedData = vtkSmartPointer < vtkAggregateDataSetFilter >::New();
  reducedData->SetInputData(output);
  reducedData->Update();
  vtkIdType numPointsBefore = reducedData->GetUnstructuredGridOutput()->GetNumberOfPoints();
  
  vtkSmartPointer < vtkCleanUnstructuredGrid > clean = vtkSmartPointer < vtkCleanUnstructuredGrid >::New(); 
  clean->SetInputData(reducedData->GetOutput());
  clean->Update();
  output->ShallowCopy(clean->GetOutput());
  vtkIdType numPointsAfter = clean->GetOutput()->GetNumberOfPoints();

  end = std::chrono::steady_clock::now();
  std::cout << "[reduceDataSet] Elapsed time in milliseconds : "
	  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
	  << " ms" << std::endl;
  std::cout << "[reduceDataSet] critical cells: " << output->GetNumberOfCells() << std::endl;
  std::cout << "[reduceDataSet] points removed: " << numPointsBefore - numPointsAfter << std::endl;

  vtkDataSetWriter * test = vtkDataSetWriter::New();
  test->SetInputData(output);
  test->SetFileName("aggregate.vtk");
  test->Update();

  start = std::chrono::steady_clock::now();
  if(output->GetNumberOfCells() > 0)
  {
	double singularitySecond[3] = {0.000, 0.000, 0.000};  
  	CriticalPointExtractor cp_cleanup(output, singularitySecond); //perturbation is ON by default
  	cp_cleanup.ComputeCriticalCells(output);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "[cleanup] Elapsed time in milliseconds : "
  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
  << " ms" << std::endl;
  std::cout << "[cleanup] Stable critical cells: " << output->GetNumberOfCells() << std::endl;

  test->SetInputData(output);
  test->SetFileName("cleaned.vtk");
  test->Update();

//   if(output->GetNumberOfCells() > 0)
//   {
// 	double singularitySecond[3] = {0.000, 0.000, 0.000};  
//   	CriticalPointExtractor cp_cleanup(output, singularitySecond, true);
//   	cp_cleanup.ComputeCriticalCells(output);
//   }
//   end = std::chrono::steady_clock::now();
//   std::cout << "[cleanup] Elapsed time in milliseconds : "
//   << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
//   << " ms" << std::endl;
//   std::cout << "[cleanup] Stable critical cells: " << output->GetNumberOfCells() << std::endl;

//   test->SetInputData(output);
//   test->SetFileName("cleaned2.vtk");
//   test->Update();
  return 1;
}

CriticalPointExtractor::CriticalPointExtractor(vtkSmartPointer<vtkDataSet> input,
											   double *currentSingularity,
											   bool pertubate)
{
	//Configure openmp
	numThreads = std::thread::hardware_concurrency() * 4; //!< Number of OpenMP threads
	omp_set_num_threads(numThreads);
	
	//Store singularity
	singularity[0] = currentSingularity[0];
	singularity[1] = currentSingularity[1];
	singularity[2] = currentSingularity[2];

	vtkIdType numPoints = input->GetNumberOfPoints();
	// ZERO_ID = 1000000000 + 1;
	ZERO_ID = numPoints + 1;

	//Allocate memory
	vecPointCoordinates.resize(numPoints);
	vecVectors.resize(numPoints);
	vecPerturbation.resize(numPoints);
	
	//Store vectors and point coordinates for internal usage
	vtkSmartPointer<vtkDataArray> vectors = input->GetPointData()->GetVectors();

	if(pertubate)
		Perturbate(singularity, ZERO_ID);

	position = new double[numPoints*3];
	vector = new double[numPoints*3];
	perturbation = new double[numPoints*3];

	#pragma omp parallel for
	for(vtkIdType i=0; i < numPoints; i++) 
	{
		vectors->GetTuple(i,&vector[i * 3]);
		input->GetPoint(i, &position[i * 3]);		
		
		vtkIdType globalID = i;//GlobalUniqueID(&position[i * 3]);
		if(pertubate)
			Perturbate(&perturbation[i * 3], globalID);

		vecVectors[i] = &vector[i * 3];
		vecPerturbation[i] = &perturbation[i * 3];
		vecPointCoordinates[i] = &position[i * 3];
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
	}

	//Allocate size for cells which depends on input cell type
	if(VTK_PIXEL == cellType || VTK_QUAD == cellType) {
		vecCellIds.resize(numCells * 2);
		numCellIds=3;
	}
	else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType) {
		vecCellIds.resize(numCells * 5);
		numCellIds=4;
	}
	else if (VTK_TRIANGLE == cellType) {
		vecCellIds.resize(numCells);
		numCellIds=3;
	}
	else if (VTK_TETRA == cellType) {
		vecCellIds.resize(numCells);
		numCellIds=4;
	}

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
			vecCellIds[i * 2]     = new vtkIdType[3]{ ids->GetId(0) , ids->GetId(1), ids->GetId(2)};
			vecCellIds[i * 2 + 1] = new vtkIdType[3]{ ids->GetId(1) , ids->GetId(3), ids->GetId(2)};
		}
		else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType)
		{
			vecCellIds[i * 5]     = new vtkIdType[4]{ ids->GetId(0) , ids->GetId(6), ids->GetId(4), ids->GetId(5)};
			vecCellIds[i * 5 + 1] = new vtkIdType[4]{ ids->GetId(3) , ids->GetId(5), ids->GetId(7), ids->GetId(6)};
			vecCellIds[i * 5 + 2] = new vtkIdType[4]{ ids->GetId(3) , ids->GetId(1), ids->GetId(5), ids->GetId(0)};
			vecCellIds[i * 5 + 3] = new vtkIdType[4]{ ids->GetId(0) , ids->GetId(3), ids->GetId(2), ids->GetId(6)};
			vecCellIds[i * 5 + 4] = new vtkIdType[4]{ ids->GetId(0) , ids->GetId(6), ids->GetId(3), ids->GetId(5)};
		}
		else if (VTK_TRIANGLE == cellType)
		{
			vecCellIds[i] = new vtkIdType[3]{ ids->GetId(0) , ids->GetId(1), ids->GetId(2)};
		}
		else if (VTK_TETRA == cellType)
		{
			vecCellIds[i] = new vtkIdType[4]{ ids->GetId(0) , ids->GetId(1), ids->GetId(2), ids->GetId(3)};
		}
		else {
			std::cout << "[CriticalPointExtractor] Error: unknown cell type " << std::endl;
				continue;
		}	
		
	}

	std::cout << "[CriticalPointExtractor::identify_critical_points] Extracted " << vecCellIds.size() << " cells" << std::endl;
}

// vtkIdType CriticalPointExtractor::GlobalUniqueID(double* pos)
// {
// 	//Function to calculate global unique id
// 	// std::hash<double> double_hash;
// 	// return (1000000*pos[0] + 10000000*pos[1] + 100000000*pos[2]);
// 	vtkIdType xD = static_cast<vtkIdType>(pos[0]);
// 	vtkIdType yD = static_cast<vtkIdType>(pos[1])<<8;
// 	vtkIdType zD = static_cast<vtkIdType>(pos[2])<<16;
// 	// std::cout<<"orig-cast: "<<xD<<" "<<yD<<" "<<zD<<std::endl;
// 	// GlobalIdType hashed = pos[0] + pos[1]*xR + pos[2]*xR*yR;
// 	vtkIdType hashed = xD + yD + zD + xD*yD + 2*yD*zD + 4*xD*zD + xD*yD*zD;
// 	return hashed;
// }

void CriticalPointExtractor::Perturbate(double* values, vtkIdType id) {
	// perturbation function f(e,i,j) = eps^2^i*delta-j
	// eps ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	//eps and delta are constant.. so I compute them one time at the beginning
	vtkIdType i = id + 1;
	double i_norm = static_cast<double>(i)/static_cast<double>(vecPerturbation.size());
	double exp_coeff = 1+i_norm*delta;

	for(int j=0; j<3; j++) {
		values[j] = std::pow(eps,std::pow(2,exp_coeff-static_cast<double>(j+1)));
	}
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::ComputeCriticalCells(vtkSmartPointer<vtkDataSet> output) 
{

	vtkSmartPointer<vtkUnstructuredGrid> outputData = vtkUnstructuredGrid::SafeDownCast(output);

	vtkIdType cells_num = vecCellIds.size();
	
	//Vector of critical cell ids
	std::vector<CriticalPoint> vecCriticalCellIDs;

	int matrixSize = 3;
	if(iExchangeIndex == 3)
		matrixSize = 4;

    //std::vector<DynamicMatrix> vecMatrices;
	std::vector<DynamicMatrix> vecMatrices;
	if(matrixSize == 4)
		vecMatrices.assign(numThreads,Eigen::Matrix4d());
	if(matrixSize == 3)
		vecMatrices.assign(numThreads,Eigen::Matrix3d());	

	//std::cout << "[CriticalPointExtractor::identify_critical_points] Matrix size(" << matrixSize << "," << matrixSize << ")"<< std::endl;
	//std::cout << "[CriticalPointExtractor::identify_critical_points] Exchange index: " << iExchangeIndex << std::endl;
	//std::cout << "[CriticalPointExtractor::identify_critical_points] Identifing critical cells "<< std::endl;

	//Check for every cell if a critical point (passed singularity as argument) exists
#pragma omp parallel 
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();	//Thread ID

		//Remove synchronization with nowait
		#pragma omp for 
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
	vtkSmartPointer<vtkDoubleArray> vectorField = vtkSmartPointer<vtkDoubleArray>::New();
	singularityType->SetName("singularityType");
	vectorField->SetName("Vec");
	vectorField->SetNumberOfComponents(3);

	//Prepare output data sequentially. Insert every critical cell to output
	for (auto cellID : vecCriticalCellIDs)
	{
		const vtkIdType* vecVertexIds = vecCellIds[cellID.id];
		vtkSmartPointer<vtkIdList> newPointIDs = vtkSmartPointer<vtkIdList>::New();	
		for (int index = 0; index < numCellIds; index++)
		{
			//Insert the type of singularity 
			singularityType->InsertNextTuple1(cellID.type);
			double pCoords[3];

			pCoords[0] = vecPointCoordinates[vecCellIds[cellID.id][index]][0];
			pCoords[1] = vecPointCoordinates[vecCellIds[cellID.id][index]][1];
			pCoords[2] = vecPointCoordinates[vecCellIds[cellID.id][index]][2];
			vtkIdType pointID = pointArray->InsertNextPoint(pCoords[0], pCoords[1], pCoords[2]);
			newPointIDs->InsertNextId(pointID);
			
			//Insert the vector field values
			vectorField->InsertNextTuple3(vecVectors[vecCellIds[cellID.id][index]][0],vecVectors[vecCellIds[cellID.id][index]][1],vecVectors[vecCellIds[cellID.id][index]][2]);
		}
		cellArray->InsertNextCell(newPointIDs);
	}

	//Add points and cells to polydata
	outputData->SetPoints(pointArray); 
	if(numCellIds == 4)
		outputData->SetCells(VTK_TETRA, cellArray);
	if(numCellIds == 3)
		outputData->SetCells(VTK_TRIANGLE, cellArray);
	outputData->GetPointData()->AddArray(singularityType);
	outputData->GetPointData()->AddArray(vectorField);
	outputData->GetPointData()->SetVectors(vectorField);

	// vtkSmartPointer < vtkCleanUnstructuredGrid > clean = vtkSmartPointer < vtkCleanUnstructuredGrid >::New(); 
    // clean->SetInputData(outputData);
    // clean->Update();
    // outputData->ShallowCopy(clean->GetOutput());
}

CriticalPointExtractor::CriticalPointType CriticalPointExtractor::PointInCell(/*const std::vector<vtkIdType> &ids*/ const vtkIdType* ids, DynamicMatrix &vecMatrix) {
	// 1. compute the initial sign of the determinant of the cell
	double targetDeterminant = ComputeDeterminant(ids, vecMatrix, false, 0);
	bool targetDirection = DeterminantCounterClockWise(targetDeterminant);
	//Check for non data values (vector is zero and determinant also) 
	if (targetDeterminant == 0)
	{
		return REGULAR;
	}
	
	double tmpDeterminant;
	bool tmpDirection;
	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0) given as i to Positive function
	for (int i = 1; i < numCellIds; i++) {
		// 2.2. compute the determinant sign again 
		tmpDeterminant   = ComputeDeterminant(ids, vecMatrix, false, i);
		tmpDirection    = DeterminantCounterClockWise(tmpDeterminant);
		
		// 2.3. check if it changes --> if so return false
		if (targetDirection != tmpDirection)
		{
			return REGULAR; // regular cell
		}
	}

	double initialDeterminant = ComputeDeterminant(ids, vecMatrix, true);	
	bool initialDirection     = DeterminantCounterClockWise(initialDeterminant);	
	
	if (initialDirection != targetDirection)
	{
		return SADDLE; // we found a saddle
	}

	return SINGULARITY; // the cell is critical, since the sign never change
}

double CriticalPointExtractor::ComputeDeterminant(
	//const std::vector<vtkIdType> &ids,
	const vtkIdType* ids,
	DynamicMatrix &vecMatrix,
	bool usePoints,
	long perturbationID
){
	
	std::array<vtkIdType, 4> tmpIds;
	int numIds = numCellIds;
	if(numIds == 3)
		tmpIds = {ids[0],ids[1],ids[2]};
	else
		tmpIds = {ids[0],ids[1],ids[2],ids[3]};

	//Exchanges every facet with the zero vector
	if (perturbationID != -1)
	{
		tmpIds[perturbationID] = ZERO_ID;
	}

	// 1. Sort and check swap operations (check)
	int swapOperations = Sort(&tmpIds[0], numIds);

	double vecValues[3];
	for (vtkIdType i = 0; i < numIds; i++) 
	{
		vtkIdType& pointID = tmpIds[i];
		if (pointID != ZERO_ID)
		{
			if(!usePoints)
			{
				vecValues[0] = vecVectors[pointID][0] + vecPerturbation[pointID][0] ;
				vecValues[1] = vecVectors[pointID][1] + vecPerturbation[pointID][1] ;
				vecValues[2] = vecVectors[pointID][2] + vecPerturbation[pointID][2] ;
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

		for (vtkIdType j = 0; j < vecMatrix.cols(); j++) {
			vecMatrix(i,j) = vecValues[j];
		}
		vecMatrix(i, iExchangeIndex) = 1;
	}

	/*
	//std::cout << " \t ######################################################################### " << std::endl;
	//std::cout << " \t\tVertex IDs ";
	//for (int x = 0; x < tmpIds.size(); x++)
	//	std::cout << tmpIds[x] << " ";
	//std::cout << std::endl;
	//std::cout << " \t\tMatrix: " << std::endl;
	//std::cout << " \t\t ";
	//for (int x = 0; x < vecMatrix.rows(); x++)
	//{
	//	for (int y = 0; y < vecMatrix.cols(); y++)
	//		std::cout << vecMatrix(x, y) << " ";
	//	std::cout << std::endl;
	//	std::cout << " \t\t ";
	//}
	//std::cout << std::endl;
	//std::cout << " \t ######################################################################### " << std::endl;
	//*/

	// 2. compute determinant sign
	double det = 0;
	if(numIds == 4)
		det = static_cast<Eigen::Matrix4d>(vecMatrix).determinant();
	else
		det = static_cast<Eigen::Matrix3d>(vecMatrix).determinant();	

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		det *= -1;
	}
	return det;
}

int CriticalPointExtractor::Sort(vtkIdType* ids, int n)
{
	if (n == 3) //Triangle
	{
		return Sort3(ids);
	}
	else if (n == 4) //TETRA 
	{
		return Sort4(ids);
	}
	else
	{
		std::cout << "Warning cell type currently not supported" << std::endl;
		return 0;
	}
}

int CriticalPointExtractor::Sort3(vtkIdType* ids)
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

int CriticalPointExtractor::Sort4(vtkIdType* ids)
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

bool CriticalPointExtractor::DeterminantCounterClockWise(double& det)
{
	if (det > 0)
		return true;
	else
		return false;
}
