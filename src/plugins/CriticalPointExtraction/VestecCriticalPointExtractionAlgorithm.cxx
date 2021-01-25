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
#include <vtkImageData.h>

#include <sstream>
#include <chrono>
#include <cmath>
#include <thread>

#include <omp.h>
#include <random>
#include <array>

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

	// Local cleanup done by every worker
	start = std::chrono::steady_clock::now();	
	//Collect all critical cells per MPI rank
	vtkSmartPointer < vtkAggregateDataSetFilter > reducedData = vtkSmartPointer < vtkAggregateDataSetFilter >::New();
	reducedData->SetInputData(output);
	reducedData->Update();
	vtkIdType numPointsBefore = reducedData->GetUnstructuredGridOutput()->GetNumberOfPoints();
	end = std::chrono::steady_clock::now();
	if(mpiRank == 0) 
	std::cout << "[MPI:" << mpiRank << "] [RequestData::reduceDataSet] Elapsed time in milliseconds : "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " ms" << std::endl;
	
	start = std::chrono::steady_clock::now();
	vtkSmartPointer < vtkCleanUnstructuredGrid > clean = vtkSmartPointer < vtkCleanUnstructuredGrid >::New(); 
	clean->SetInputData(reducedData->GetOutput());
	clean->Update();
	output->ShallowCopy(clean->GetOutput());
	vtkIdType numPointsAfter = clean->GetOutput()->GetNumberOfPoints();

	end = std::chrono::steady_clock::now();
	// if(mpiRank == 0) 
	if(output->GetNumberOfCells() > 0)
	{
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] Elapsed time in milliseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] critical cells: " << output->GetNumberOfCells() << std::endl;
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] points removed: " << numPointsBefore - numPointsAfter << std::endl;
	}

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

CriticalPointExtractor::CriticalPointExtractor(vtkSmartPointer<vtkDataSet> input,
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

	//Allocate memory
	vecPointCoordinates.resize(numPoints);
	vecVectors.resize(numPoints);
	vecPerturbation.resize(numPoints);
	
	//Store vectors and point coordinates for internal usage
	vtkSmartPointer<vtkDataArray> vectors = input->GetPointData()->GetVectors();	

	position = new double[numPoints*3];
	vector = new double[numPoints*3];
	perturbation = new double[numPoints*3];
	
	//Get the local bounds of the current MPI process
	double local_bounds[6];
	input->GetBounds(local_bounds);
	double xDim = fabs(local_bounds[1] - local_bounds[0]);
	double yDim = fabs(local_bounds[3] - local_bounds[2]);
	double zDim = fabs(local_bounds[5] - local_bounds[4]);	

	/// then extract some global characteristics of the dataset
	/// like, 1. global bounds
	double global_bounds[6];
	vtkSmartPointer<vtkMultiProcessController> controller = vtkMultiProcessController::GetGlobalController();
	controller->AllReduce(&local_bounds[0], &global_bounds[0], 1, vtkCommunicator::StandardOperations::MIN_OP);
	controller->AllReduce(&local_bounds[2], &global_bounds[2], 1, vtkCommunicator::StandardOperations::MIN_OP);
	controller->AllReduce(&local_bounds[4], &global_bounds[4], 1, vtkCommunicator::StandardOperations::MIN_OP);
	controller->AllReduce(&local_bounds[1], &global_bounds[1], 1, vtkCommunicator::StandardOperations::MAX_OP);
	controller->AllReduce(&local_bounds[3], &global_bounds[3], 1, vtkCommunicator::StandardOperations::MAX_OP);
	controller->AllReduce(&local_bounds[5], &global_bounds[5], 1, vtkCommunicator::StandardOperations::MAX_OP);

	int mpiRanks = controller->GetNumberOfProcesses(); // to check how many MPI processes are up and running

	// 2. global sides
	double global_sides[3] = { fabs(global_bounds[1]-global_bounds[0]), fabs(global_bounds[3]-global_bounds[2]), fabs(global_bounds[5]-global_bounds[4]) };
	// 3. global maximum coordinates
	double global_max_coords[3] = { global_bounds[1] , global_bounds[3] , global_bounds[5] };
	// 4. spacing between points (WARNING: this works only on regular grids!!)	
	double spacing[3];
	vtkImageData::SafeDownCast(input)->GetSpacing(spacing);
	// 5. and the global extent (i.e., the resolution of the dataset)
	int global_extent[6] = { 0, global_sides[0]/spacing[0], 0, global_sides[1]/spacing[1], 0, global_sides[2]/spacing[2]};

	// the global max id is needed for computing the perturbation in each point
	// -- if we have just one MPI process then we can directly use the number of points value, since the indexing is given and consistent
	// -- otherwise, in case of multiple MPI processes we have to derive the global id from some geometric information linked to the grid
	long max_global_id = mpiRanks == 1 ? numPoints : GlobalUniqueID(global_max_coords,spacing,global_extent,global_bounds);	

	ZERO_ID = max_global_id + 1; // this is the of the singularity vector
	/// NOTICE: we do not have to perturbate the singularity vector
	if(pertubate) {
		double zero_vec[3] = {0,0,0};
	  	// Perturbate(sing_pert, ZERO_ID, max_global_id);
		// toFixed(sing_pert);

		Perturbate(zero_vec, ZERO_ID, max_global_id);		
		singularity[0] += zero_vec[0];
		singularity[1] += zero_vec[1];
		singularity[2] += zero_vec[2];
	 	// toFixed(singularity);	

		//  int a; cin>>a;	 
	}

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

	#pragma omp parallel firstprivate(vecCellPerThread)
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();//Thread ID

		#pragma omp for nowait
		for(vtkIdType i=0; i < numPoints; i++) 
		{
			vectors->GetTuple(i,&vector[i * 3]);
			input->GetPoint(i, &position[i * 3]);	

			// -- if we have just one MPI process then we can directly use the point id, since the indexing is given and consistent
			// -- otherwise, in case of multiple MPI processes we have to derive the global id of the point from some geometric information linked to the grid
			long global_id = mpiRanks == 1 ? i : GlobalUniqueID(&position[i * 3],spacing,global_extent,global_bounds);

			/// ADD CNL CALLS HERE
			// toFixed(&vector[i*3]);
			// toFixed(&position[i*3]);		
			
			if(pertubate) {								
				Perturbate(&perturbation[i * 3], global_id, max_global_id);				
				// toFixed(&perturbation[i*3]);
			}
			
			// int a; cin>>a;

			vecVectors[i] 		= &vector[i * 3];
			vecPerturbation[i] 	= &perturbation[i * 3];
			vecPointCoordinates[i] 	= &position[i * 3];
		}

		#pragma omp for
		for (vtkIdType i = 0; i < numCells; i++) {
			//get the current cell
			input->GetCell(i, vecCellPerThread[threadIdx]);

			//Get the associated point ids for the cell 
			vtkIdList* ids = vecCellPerThread[threadIdx]->GetPointIds();
			
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
				std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Error: unknown cell type " << std::endl;
					continue;
			}	
			
		}
	}

	#pragma omp parallel for
    for(int x=0; x < numThreads;++x)
	{
		vecCellPerThread[x]->Delete();
	}
	vecCellPerThread.clear();
	if(mpiRank == 0) std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Extracted " << vecCellIds.size() << " simplices" << std::endl;
}

long CriticalPointExtractor::GlobalUniqueID(double* pos, double *spacing, int *global_extent, double * global_bounds)
{
	///Function that calculates global unique id

	/// 1. structured coordinates
	long x = std::lround(pos[0]/spacing[0]-global_bounds[0]);
	long y = std::lround(pos[1]/spacing[1]-global_bounds[2]);
	long z = std::lround(pos[2]/spacing[2]-global_bounds[4]);

	/// 2. then compute the resolution
	long resx = global_extent[1]+1;
	long resy = global_extent[3]+1;
	long resz = global_extent[5]+1;

	/// 3. then the global id
	long globalid = z * resy * resx + y * resz + x;

	// x + y * xDim + z * xDim * yDim
	// return fabs(pos[0]) + fabs(pos[1]) * global_dims[0] + fabs(pos[2]) * global_dims[0] * global_dims[1];
	return globalid;
}

void CriticalPointExtractor::toFixed(double *values) {	
	// std::cout<<"original: "<<values[0]<<" "<<values[1]<<" "<<values[2]<<std::endl;
	
	auto tofixed = cnl::scaled_integer<long long, cnl::power<-14>>(values[0]);
	values[0] = to_rep(tofixed);

	tofixed = cnl::scaled_integer<long long, cnl::power<-14>>(values[1]);
	values[1] = to_rep(tofixed);

	tofixed = cnl::scaled_integer<long long, cnl::power<-14>>(values[2]);
	values[2] = to_rep(tofixed);

	// std::cout<<"to_fixed: "<<values[0]<<" "<<values[1]<<" "<<values[2]<<std::endl;
}

void CriticalPointExtractor::Perturbate(double* values, long id, long max_global_id) {
	// perturbation function f(e,i,j) = eps^2^i*delta-j
	// eps ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	//eps and delta are constant.. so I compute them one time at the beginning
	// vtkIdType i = id + 1;
	double i_norm = static_cast<double>(id)/static_cast<double>(max_global_id);
	double exp_coeff = 1+i_norm*delta;

	for(int j=0; j<3; j++) {
		values[j] = std::pow(eps,std::pow(2,exp_coeff-static_cast<double>(j+1)));
	}
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::ComputeCriticalCells() 
{
	vtkIdType cells_num = vecCellIds.size();
	
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
#pragma omp parallel firstprivate(vecMatrices)
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
}

void CriticalPointExtractor::writeCriticalCells(vtkSmartPointer<vtkDataSet> output) 
{
	vtkSmartPointer<vtkUnstructuredGrid> outputData = vtkUnstructuredGrid::SafeDownCast(output);

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
