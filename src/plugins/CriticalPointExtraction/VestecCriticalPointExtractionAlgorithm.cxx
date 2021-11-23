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
#include <vtkMergePoints.h>
#include <vtkDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkAppendFilter.h>

#include <sstream>
#include <chrono>
#include <cmath>
#include <thread>

#include <omp.h>
#include <random>
#include <array>
#include <functional> //for std::hash

#include <Eigen/Core>


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

	vtkSmartPointer<vtkMultiProcessController> controller = vtkMultiProcessController::GetGlobalController();
	int mpiRank  = controller->GetLocalProcessId();
	int mpiRanks = controller->GetNumberOfProcesses();

	//Set active array
	vtkDataArray *inScalars = this->GetInputArrayToProcess(0, inputVector);
	input->GetPointData()->SetActiveVectors(inScalars->GetName());
	
	//Compute critical points
	double singularity[3] = {0.000, 0.000, 0.000}; 
	
	controller->Barrier();

	auto start = std::chrono::steady_clock::now();
	CriticalPointExtractor cp_extractor(input, singularity, mpiRank); //perturbation is ON by default
	auto end = std::chrono::steady_clock::now();
	if(mpiRank == 0)  	
	std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::Constructor] Elapsed time in milliseconds : "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " ms" << std::endl;
	double constructor_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	controller->Barrier();
	
	start = std::chrono::steady_clock::now();
	cp_extractor.ComputeCriticalCells();
	end = std::chrono::steady_clock::now();
	if(mpiRank == 0) 
	std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Elapsed time in milliseconds : "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " ms" << std::endl;
	if(mpiRank == 0) 
	std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Unstable critical cells: " << output->GetNumberOfCells() << std::endl;
	double critical_point_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

	//Write cells to output
	cp_extractor.writeCriticalCells(output);	

	//Sync to measure correct time for the reduce
	controller->Barrier();

	int receiveProc = 0;
  	std::vector<vtkSmartPointer<vtkDataObject>> recvBuffer;
	start = std::chrono::steady_clock::now();	
  	controller->Gather(output, recvBuffer, receiveProc);
	end = std::chrono::steady_clock::now();
	if(mpiRank == receiveProc) 
	{
		std::cout << "[MPI:" << mpiRank << "] [RequestData::reduceDataSet] Elapsed time in milliseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;
	}	
	
	//Sync to measure correct time for the clean step
	controller->Barrier();
	
	start = std::chrono::steady_clock::now();	
	if (mpiRank == receiveProc)
	{
		if (recvBuffer.size() != 0)
		{
			vtkNew<vtkAppendFilter> appendFilter;
			appendFilter->MergePointsOn();
			for (std::vector<vtkSmartPointer<vtkDataObject> >::iterator it = recvBuffer.begin();
				it != recvBuffer.end(); ++it)
			{
				appendFilter->AddInputData(*it);
			}
			appendFilter->Update();
			output->ShallowCopy(appendFilter->GetOutput());
		}
 	}	
	end = std::chrono::steady_clock::now();

	//Collect all degenerated cases
	controller->AllReduce(&cp_extractor.local_deg_cases, &cp_extractor.global_deg_cases, 1, vtkCommunicator::StandardOperations::SUM_OP);
	if(mpiRank == receiveProc) 
	{
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] Elapsed time in milliseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] critical cells: " << output->GetNumberOfCells() << std::endl;
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] degenerate cases: " << cp_extractor.global_deg_cases << std::endl;
	}

	//For benchmarking we aggregate all timinigs from the mpi processes
	if(mpiRanks > 1) 
	{
		// aggregating timings
		double tot_constructor_time=0;
		controller->Reduce(&constructor_time, &tot_constructor_time, 1, vtkCommunicator::StandardOperations::MAX_OP, 0);
		if(mpiRank == 0) 
			std::cout<<"[MPI:" << mpiRank << "] [RequestData] [MAX] Constructor time: "<<tot_constructor_time<<std::endl;
		double tot_cp_time=0;
		controller->Reduce(&critical_point_time,&tot_cp_time,1,vtkCommunicator::StandardOperations::MAX_OP,0);
		if(mpiRank == 0) 
			std::cout<<"[MPI:" << mpiRank << "] [RequestData] [MAX] Critical Point Extraction time: "<<tot_cp_time<<std::endl;		
	}
	
	return 1;
}

CriticalPointExtractor::CriticalPointExtractor(vtkDataSet* input,
											   double *currentSingularity, int mpiRank,
											   bool pertubate)
{
	//Configure openmp
	numThreads = omp_get_max_threads(); //!< Number of OpenMP threads
	if(mpiRank == 0) std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor] Number of threads: " << numThreads << std::endl;
	
	//Store singularity
	singularity[0] = currentSingularity[0];
	singularity[1] = currentSingularity[1];
	singularity[2] = currentSingularity[2];

	vtkIdType numPoints = input->GetNumberOfPoints();	

	//Store vectors and point coordinates for internal usage
	vtkDataArray* vectors = input->GetPointData()->GetVectors();	

	position = new double[numPoints*3];
	vector = new double[numPoints*3];

	//Data that needs to stored in memory
	//Cell indices, vertex positions and velocity vector
	//Memory for cell indices is added later
	long long max_memory = (numPoints*6)*sizeof(double);
	
	//Store metadata such as bounds, extents, dimension, ...
	DataSetMetadata dm;	
	
	//Get the local bounds of the current MPI process
	input->GetBounds(dm.local_bounds);			

	// get global lenght of each axis
	double xDim = fabs(dm.local_bounds[1] - dm.local_bounds[0]);
	double yDim = fabs(dm.local_bounds[3] - dm.local_bounds[2]);
	double zDim = fabs(dm.local_bounds[5] - dm.local_bounds[4]);	

	// this is the of the singularity vector
	ZERO_ID = dm.max_global_id; 

	iExchangeIndex = 3; 					//3D dataset 
	if (xDim == 0.0) iExchangeIndex = 0; 	//2D dataset with yz
	if (yDim == 0.0) iExchangeIndex = 1; 	//2D dataset with xz
	if (zDim == 0.0) iExchangeIndex = 2; 	//2D dataset with xy

	//Configure for parallel independent processing
    std::vector<vtkGenericCell*> vecCellPerThread;      //Cell for each thread

	//Get the cell type once (needed for correct allocation)
	vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
	vtkIdType numCells = input->GetNumberOfCells();
	input->GetCell(0, cell);
    vtkIdType cellType = cell->GetCellType();
	
	vecCellPerThread.resize(numThreads);
	#pragma omp parallel for
    for(int x=0; x < numThreads;++x)
	{
		vecCellPerThread[x] = vtkGenericCell::New();
	}

	//Allocate size for cells which depends on input cell type
	if(VTK_PIXEL == cellType || VTK_QUAD == cellType) {
		vecCellIds = new vtkIdType[numCells * 6];
		numCellIds=3;
		numSimplices = numCells * 2;
		numSimplicesPerCell = 2;
	}
	else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType) {
		vecCellIds = new vtkIdType[numCells * 20];
		numCellIds=4;
		numSimplices = numCells * 5;
		numSimplicesPerCell = 5;
	}
	else if (VTK_TRIANGLE == cellType) {
		vecCellIds = new vtkIdType[numCells * 3];
		numCellIds=3;
		numSimplices = numCells;
		numSimplicesPerCell = 1;
	}
	else if (VTK_TETRA == cellType) {
		vecCellIds = new vtkIdType[numCells * 4];
		numCellIds=4;
		numSimplices = numCells;
		numSimplicesPerCell = 1;
	}	

	#pragma omp parallel firstprivate(vecCellPerThread) //private(dm,mpiRanks)
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();//Thread ID

		#pragma omp for nowait
		for(vtkIdType i=0; i < numPoints; i++) 
		{
			vectors->GetTuple(i,&vector[i * 3]);
			input->GetPoint(i, &position[i * 3]);

			// We use a unique hash for perturbating the vector field
			vtkIdType hash = ComputeHash(&position[i * 3]);
			Perturbate(&vector[i * 3], hash, dm);
		}
	
		#pragma omp for
		for (vtkIdType i = 0; i < numCells; i++) {
			//get the current cell
			input->GetCell(i, vecCellPerThread[threadIdx]); /// input is a vtkSmartPointer

			//Get the associated point ids for the cell 
			vtkIdList* ids = vecCellPerThread[threadIdx]->GetPointIds();
			
			if (VTK_PIXEL == cellType || VTK_QUAD == cellType)
			{
				vecCellIds[i*6] = ids->GetId(0); vecCellIds[i*6+1] = ids->GetId(1); vecCellIds[i*6+2] = ids->GetId(2);
				vecCellIds[i*6+3] = ids->GetId(1); vecCellIds[i*6+4] = ids->GetId(3); vecCellIds[i*6+5] = ids->GetId(2);
			}
			else if (VTK_VOXEL == cellType)
			{
				vecCellIds[i*20] = ids->GetId(0); vecCellIds[i*20+1] = ids->GetId(6); vecCellIds[i*20+2] = ids->GetId(4); vecCellIds[i*20+3] = ids->GetId(5);
				vecCellIds[i*20+4] = ids->GetId(3); vecCellIds[i*20+5] = ids->GetId(5); vecCellIds[i*20+6] = ids->GetId(7); vecCellIds[i*20+7] = ids->GetId(6);
				vecCellIds[i*20+8] = ids->GetId(3); vecCellIds[i*20+9] = ids->GetId(1); vecCellIds[i*20+10] = ids->GetId(5); vecCellIds[i*20+11] = ids->GetId(0);
				vecCellIds[i*20+12] = ids->GetId(0); vecCellIds[i*20+13] = ids->GetId(3); vecCellIds[i*20+14] = ids->GetId(2); vecCellIds[i*20+15] = ids->GetId(6);
				vecCellIds[i*20+16] = ids->GetId(0); vecCellIds[i*20+17] = ids->GetId(6); vecCellIds[i*20+18] = ids->GetId(3); vecCellIds[i*20+19] = ids->GetId(5);
			}
			else if (VTK_HEXAHEDRON == cellType)
			{
				if(i % 2){
					vecCellIds[i*20] = ids->GetId(0); vecCellIds[i*20+1] = ids->GetId(1); vecCellIds[i*20+2] = ids->GetId(3); vecCellIds[i*20+3] = ids->GetId(4);
					vecCellIds[i*20+4] = ids->GetId(1); vecCellIds[i*20+5] = ids->GetId(4); vecCellIds[i*20+6] = ids->GetId(5); vecCellIds[i*20+7] = ids->GetId(6);
					vecCellIds[i*20+8] = ids->GetId(1); vecCellIds[i*20+9] = ids->GetId(4); vecCellIds[i*20+10] = ids->GetId(6); vecCellIds[i*20+11] = ids->GetId(3);
					vecCellIds[i*20+12] = ids->GetId(1); vecCellIds[i*20+13] = ids->GetId(3); vecCellIds[i*20+14] = ids->GetId(6); vecCellIds[i*20+15] = ids->GetId(2);
					vecCellIds[i*20+16] = ids->GetId(3); vecCellIds[i*20+17] = ids->GetId(6); vecCellIds[i*20+18] = ids->GetId(7); vecCellIds[i*20+19] = ids->GetId(4);
				} else {
					vecCellIds[i*20] = ids->GetId(2); vecCellIds[i*20+1] = ids->GetId(1); vecCellIds[i*20+2] = ids->GetId(5); vecCellIds[i*20+3] = ids->GetId(0);
					vecCellIds[i*20+4] = ids->GetId(0); vecCellIds[i*20+5] = ids->GetId(2); vecCellIds[i*20+6] = ids->GetId(3); vecCellIds[i*20+7] = ids->GetId(7);
					vecCellIds[i*20+8] = ids->GetId(2); vecCellIds[i*20+9] = ids->GetId(5); vecCellIds[i*20+10] = ids->GetId(6); vecCellIds[i*20+11] = ids->GetId(7);
					vecCellIds[i*20+12] = ids->GetId(0); vecCellIds[i*20+13] = ids->GetId(7); vecCellIds[i*20+14] = ids->GetId(4); vecCellIds[i*20+15] = ids->GetId(5);
					vecCellIds[i*20+16] = ids->GetId(0); vecCellIds[i*20+17] = ids->GetId(2); vecCellIds[i*20+18] = ids->GetId(7); vecCellIds[i*20+19] = ids->GetId(5);
				}
			}
			else if (VTK_TRIANGLE == cellType)
			{
				vecCellIds[i*3] = ids->GetId(0); vecCellIds[i*3+1] = ids->GetId(1); vecCellIds[i*3+2] = ids->GetId(2);
			}
			else if (VTK_TETRA == cellType)
			{
				vecCellIds[i*4] = ids->GetId(0); vecCellIds[i*4+1] = ids->GetId(1); vecCellIds[i*4+2] = ids->GetId(2); vecCellIds[i*4+3] = ids->GetId(3);
			}
			else {
				std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Error: unknown cell type! Type is "<< cellType << std::endl;
				continue;
			}	
			
		}
	}

	max_memory += numSimplices*numCellIds*sizeof(vtkIdType);
	if(mpiRank == 0) std::cout<<"Allocated memory: "<< max_memory / std::pow(1024,2) <<" (MB) " << std::endl;	

	// /// ==== DEBUG ONLY === ///
	// /// write to 2 separate files the points and simplexes arrays
	// std::ofstream points_file;
	// points_file.open("points.txt");
	// for(vtkIdType i=0; i < numPoints; i++) 
	// {
	// 	points_file << position[i*3] << " " << position[i*3+1] << " " << position[i*3+2] << " ";
	// 	points_file << vector[i*3] << " " << vector[i*3+1] << " " << vector[i*3+2] << std::endl;
	// }
	// points_file.close();
	// std::ofstream simplices_file;
	// simplices_file.open("simplices.txt");
	// for(vtkIdType i=0; i < vecCellIds.size(); i++) 
	// {
	// 	simplices_file << vecCellIds[i][0] << " " << vecCellIds[i][1] << " " << vecCellIds[i][2] << " " << vecCellIds[i][3] << std::endl;		
	// }
	// simplices_file.close();
	// /// ==== DEBUG ONLY === ///

	#pragma omp parallel for
    for(int x=0; x < numThreads;++x)
	{
		vecCellPerThread[x]->Delete();
	}
	vecCellPerThread.clear();
	
	std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Extracted " << numSimplices << " simplices" << std::endl;
}

vtkIdType CriticalPointExtractor::ComputeHash(double* pos) 
{
	vtkIdType globalid;
	std::hash<double*> constexpr h;
	globalid = h(pos);
	return globalid;
}

void CriticalPointExtractor::Perturbate(double* values, vtkIdType &id, DataSetMetadata &dm) {
	// perturbation function f(e,i,j) = eps^2^i*delta-j
	// eps ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	//eps and delta are constant.. so I compute them one time at the beginning
	vtkIdType i = id + 1;
	double i_norm = (static_cast<double>(i)-static_cast<double>(dm.min_global_id))/(static_cast<double>(dm.max_global_id)-static_cast<double>(dm.min_global_id));
	double exp_coeff = i_norm*delta;

	double j_norm;

	for(int j=0; j<3; j++) {
		j_norm = static_cast<double>(j+1)/3; //since we are normalizing the point id, we need to normalize as well the j-id --> to keep the perturbation small
		values[j] += std::pow(eps,std::pow(2,exp_coeff-j_norm));
	}
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::ComputeCriticalCells() 
{
	int matrixSize = 3;
	if(iExchangeIndex == 3)
		matrixSize = 4;

    std::vector<DynamicMatrix> vecMatrices;
    vecMatrices.resize(numThreads);
	if(matrixSize == 4)
		vecMatrices.assign(numThreads,Eigen::Matrix4d());
	if(matrixSize == 3)
		vecMatrices.assign(numThreads,Eigen::Matrix3d());	

	//Check for every cell if a critical point (passed singularity as argument) exists
	#pragma omp parallel firstprivate(vecMatrices) //private(vector,position)//private(singularity,ZERO_ID)
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();	//Thread ID

		//Remove synchronization with nowait
		#pragma omp for
		for (vtkIdType i = 0; i < numSimplices; i++) {
			//If the cell contains a the singularity add them to the output and we can break
			PointType ret = PointInCell(&vecCellIds[i*numCellIds], vecMatrices[threadIdx]);
			if (ret != REGULAR_POINT) {
				/// classify the critical simplex by computing the eigenvalues on its Jacobian
				ret = ClassifyCriticalSimplex(&vecCellIds[i*numCellIds]);
				CriticalPoint tmp(i,ret);
#pragma omp critical
				vecCriticalCellIDs.push_back(tmp);
			}
		}
	}
}

void CriticalPointExtractor::writeCriticalCells(vtkSmartPointer<vtkDataSet> output) 
{
	vtkSmartPointer<vtkUnstructuredGrid> outputData = vtkUnstructuredGrid::SafeDownCast(output);

	vtkSmartPointer<vtkPoints> pointArray = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIntArray> singularityType = vtkSmartPointer<vtkIntArray>::New();
	singularityType->SetName("singularityType");

	//Prepare output data sequentially. Insert every critical cell to output
	for (auto cellID : vecCriticalCellIDs)
	{
		vtkIdType positionInArray = cellID.id * numCellIds;
		const vtkIdType* vecVertexIds = &vecCellIds[positionInArray];
		vtkSmartPointer<vtkIdList> newPointIDs = vtkSmartPointer<vtkIdList>::New();	

		double* barycenter = ComputeBarycentricCoordinates(vecVertexIds);
		singularityType->InsertNextTuple1(cellID.type);

		vtkIdType newPointID = pointArray->InsertNextPoint(barycenter[0], barycenter[1], barycenter[2]);
		newPointIDs->InsertNextId(newPointID);
		cellArray->InsertNextCell(newPointIDs);
	}
		
	//Add points and cells to polydata
	outputData->SetPoints(pointArray); 
	outputData->SetCells(VTK_VERTEX, cellArray);
	outputData->GetPointData()->AddArray(singularityType);
}

CriticalPointExtractor::PointType CriticalPointExtractor::PointInCell(const vtkIdType* ids, DynamicMatrix &vecMatrix) {		
	std::array<vtkIdType, 4> tmpIds;
	int numIds = numCellIds;
	if(numIds == 3)
		tmpIds = {ids[0],ids[1],ids[2]};
	else
		tmpIds = {ids[0],ids[1],ids[2],ids[3]};

	// 1. compute the initial sign of the determinant of the cell
	double targetDeterminant = ComputeDeterminant(tmpIds, vecMatrix, false, 0);
	bool targetDirection = DeterminantCounterClockWise(targetDeterminant);
	//Check for non data values (vector is zero and determinant also) 
	if (targetDeterminant == 0)
	{
		#pragma omp atomic		
		this->local_deg_cases++;
		return REGULAR_POINT;
	}
	
	double tmpDeterminant;
	bool tmpDirection;
	// 2. for each facet (i.e. an edge in a triangle or a triangle in a tetrahedron) do
	// 2.1. replace each row of the matrix with the origin vector (0,0) or (0,0,0) given as i to Positive function
	for (int i = 1; i < numCellIds; i++) {
		if(numIds == 3)
			tmpIds = {ids[0],ids[1],ids[2]};
		else
			tmpIds = {ids[0],ids[1],ids[2],ids[3]};
		// 2.2. compute the determinant sign again 
		tmpDeterminant   = ComputeDeterminant(tmpIds, vecMatrix, false, i);
		tmpDirection    = DeterminantCounterClockWise(tmpDeterminant);
		
		// 2.3. check if it changes --> if so return false
		if (targetDirection != tmpDirection)
		{
			return REGULAR_POINT; // regular cell
		}
	}
	return UNCLASSIFIED_SINGULARITY; // the cell is critical, since the sign never change
}

CriticalPointExtractor::PointType CriticalPointExtractor::ClassifyCriticalSimplex(const vtkIdType* ids/*, DynamicMatrix &jacobian*/) {
	int numIds = numCellIds;

	CriticalPointExtractor::PointType ret = UNCLASSIFIED_SINGULARITY;
	int posReal = 0;
	int negReal = 0;

	int zeroImag = 0;
	int complexImag = 0;

	if(numIds == 4) {//in 3D for tetrahedral meshes we have a 3x3 matrix,
		Eigen::Matrix3d coordsMatrix;
		Eigen::Matrix3d vectorsMatrix;
		InitializeMatrices(ids,coordsMatrix,vectorsMatrix);

		Eigen::Matrix3d jacobianMatrix;
		jacobianMatrix = coordsMatrix * vectorsMatrix.inverse();
		auto ev = jacobianMatrix.inverse().eigenvalues();
		CheckEigenvalues(ev,posReal,negReal,zeroImag,complexImag);
		
		ret = GetCriticalSimplexType(ev,posReal,negReal,zeroImag,complexImag);		
	} else {//while in 2D for triangle meshes we have a 2x2 matrix
		Eigen::Matrix2d coordsMatrix;
		Eigen::Matrix2d vectorsMatrix;
		InitializeMatrices(ids,coordsMatrix,vectorsMatrix);
		
		Eigen::Matrix2d jacobianMatrix;
		jacobianMatrix = coordsMatrix * vectorsMatrix.inverse();
		auto ev = jacobianMatrix.inverse().eigenvalues();
		CheckEigenvalues(ev,posReal,negReal,zeroImag,complexImag);

		ret = GetCriticalSimplexType(ev,posReal,negReal,zeroImag,complexImag);
	}

	return ret;
}

double* CriticalPointExtractor::ComputeBarycentricCoordinates(const vtkIdType* ids) {
	// int numIds = numCellIds;	
	vtkIdType numDims = (numCellIds == 4) ? 3 : 2;	
	
	Eigen::MatrixXd vectorsMatrix(numDims,numDims);
	double* baricenter = new double[3]{0,0,0};
	double* lambda = new double[numDims];

	// here we initialize the vector matrix as v[i] - v[numDims]
  	for (int i = 0; i < numDims; i++) {		
		for (int j = 0; j < numDims; j++)
			vectorsMatrix(i,j) = vector[ids[j]*3+i] - vector[ids[numDims]*3+i];
	}
	vectorsMatrix = vectorsMatrix.inverse();
	// here we compute the barycentric coordinates on zero --> lambda] = T^-1(i) * vector[id[numDims]]
	for (int i = 0; i < numDims; i++) {
		lambda[i] = 0;
		for (int j = 0; j < numDims; j++)
			lambda[i] += vectorsMatrix(i,j) * vector[ids[numDims]*3+j];
	}	
	// computed the multiplication coefficient for the last vertex of the triangle/tetrahedron
	double mult_coeff = 1;
	for (int i = 0; i < numDims; i++)
		mult_coeff -= lambda[i];
	// now we have to compute the barycentric interpolation to get the coordinates
	// (in 2D) b[i] = lambda[0]*pos[0][i] + lambda[1]*pos[1][i] + pos[2][i] * mult_coeff
	// (in 3D) b[i] = lambda[0]*pos[0][i] + lambda[1]*pos[1][i] + lambda[2]*pos[2][i] + pos[3][i] * mult_coeff
	for (int i = 0; i < numDims; i++) {
		for (int j = 0; j < numDims; j++)
			baricenter[i] += position[ids[j]*3+i] * lambda[j];
		baricenter[i] += position[ids[numDims]*3+i] * mult_coeff;
	}

	return baricenter;
}

double CriticalPointExtractor::ComputeDeterminant(	
	std::array<vtkIdType, 4> &tmpIds,
	DynamicMatrix &vecMatrix,
	bool usePoints,
	int perturbationID	
){
		
	int numIds = numCellIds;
	
	//Exchanges every facet with the zero vector
	if (perturbationID != -1)
	{
		tmpIds[perturbationID] = ZERO_ID;
	}

	// 1. Sort and check swap operations (check)
	int swapOperations = Sort(&tmpIds[0], numIds);

	// double vecValues[3];
	for (vtkIdType i = 0; i < numIds; i++) 
	{
		vtkIdType& pointID = tmpIds[i];
		if (pointID != ZERO_ID)
		{
			if(!usePoints)
			{
				vecMatrix(i,0) = vector[pointID*3];
				vecMatrix(i,1) = vector[pointID*3+1];
				vecMatrix(i,2) = vector[pointID*3+2];;
			}
			else
			{
				vecMatrix(i,0) = position[pointID*3];
				vecMatrix(i,1) = position[pointID*3+1];
				vecMatrix(i,2) = position[pointID*3+2];
			}
		}
		else
		{
			vecMatrix(i,0) = singularity[0];
			vecMatrix(i,1) = singularity[1];
			vecMatrix(i,2) = singularity[2];
		}
		vecMatrix(i, iExchangeIndex) = 1;
	}

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
