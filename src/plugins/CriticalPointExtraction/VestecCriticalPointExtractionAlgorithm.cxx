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
	mpf_set_default_prec(256);	

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

	// Local cleanup done by every worker
	
	
	//Collect all critical cells per MPI rank
	//vtkSmartPointer<vtkUnstructuredGrid> gatheredOutput = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//std::vector<vtkIdType> pointCount(mpiRanks, 0);
  	//vtkIdType numPoints = output->GetNumberOfPoints();
  	//controller->AllGather(&numPoints, &pointCount[0], 1);

	int receiveProc = 0;
	/*vtkIdType maxVal = 0;
	for (int i = 0; i < mpiRanks; i++)
	{
		if (pointCount[i] > maxVal)
		{
		maxVal = pointCount[i];
		receiveProc = i;
		}
	}*/

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
		if (recvBuffer.size() > 1 && output->IsA("vtkUnstructuredGrid"))
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
	
	if(mpiRank == receiveProc) 
	{
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] Elapsed time in milliseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] critical cells: " << output->GetNumberOfCells() << std::endl;
		std::cout << "[MPI:" << mpiRank << "] [RequestData::cleanupDataSet] degenerate cases: " << cp_extractor.deg_cases << std::endl;
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

CriticalPointExtractor::CriticalPointExtractor(vtkDataSet* input,
											   double *currentSingularity, int mpiRank/*,
											   bool pertubate*/)
{
	eps = 1 / std::pow(10,14);

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

	position = new mpf_float_custom[numPoints*3];
	vector = new mpf_float_custom[numPoints*3];

	//Data that needs to be written
	long long max_memory = (numPoints*3*3)*sizeof(double);
	
	DataSetMetadata dm;	
	
	//Get the local bounds of the current MPI process
	// double local_bounds[6];
	input->GetBounds(dm.local_bounds);			

	/// then extract some global characteristics of the dataset
	/// like, 1. global bounds
	// double global_bounds[6];
	vtkSmartPointer<vtkMultiProcessController> controller = vtkMultiProcessController::GetGlobalController();
	controller->AllReduce(&dm.local_bounds[0], &dm.global_bounds[0], 1, vtkCommunicator::StandardOperations::MIN_OP);
	controller->AllReduce(&dm.local_bounds[2], &dm.global_bounds[2], 1, vtkCommunicator::StandardOperations::MIN_OP);
	controller->AllReduce(&dm.local_bounds[4], &dm.global_bounds[4], 1, vtkCommunicator::StandardOperations::MIN_OP);
	controller->AllReduce(&dm.local_bounds[1], &dm.global_bounds[1], 1, vtkCommunicator::StandardOperations::MAX_OP);
	controller->AllReduce(&dm.local_bounds[3], &dm.global_bounds[3], 1, vtkCommunicator::StandardOperations::MAX_OP);
	controller->AllReduce(&dm.local_bounds[5], &dm.global_bounds[5], 1, vtkCommunicator::StandardOperations::MAX_OP);

	int mpiRanks = controller->GetNumberOfProcesses(); // to check how many MPI processes are up and running

	// get global sides
	double xDim = fabs(dm.global_bounds[1] - dm.global_bounds[0]);
	double yDim = fabs(dm.global_bounds[3] - dm.global_bounds[2]);
	double zDim = fabs(dm.global_bounds[5] - dm.global_bounds[4]);	
	// std::cout << "xDim " << xDim << " yDim " << yDim << " zDim " << zDim << std::endl;

	// 2. global sides
	// double global_sides[3] = { xDim, yDim, zDim };
	// 3. global maximum coordinates
	mpf_float_custom global_max_coords[3] = { dm.global_bounds[1] , dm.global_bounds[3] , dm.global_bounds[5] };
	// 4. spacing between points (WARNING: this works only on regular grids!!)	
	// double spacing[3];
	if(vtkImageData::SafeDownCast(input))
		vtkImageData::SafeDownCast(input)->GetSpacing(dm.spacing);
	// 5. and the global extent (i.e., the resolution of the dataset)
	dm.global_extent = new int[6]{0, static_cast<int>(xDim/dm.spacing[0]), 0, static_cast<int>(yDim/dm.spacing[1]), 0, static_cast<int>(zDim/dm.spacing[2])};

	// int dimensions[3];
	if(vtkImageData::SafeDownCast(input))
		vtkImageData::SafeDownCast(input)->GetDimensions(dm.dimensions);

	// the global max id is needed for computing the perturbation in each point
	// -- if we have just one MPI process then we can directly use the number of points value, since the indexing is given and consistent
	// -- otherwise, in case of multiple MPI processes we have to derive the global id from some geometric information linked to the grid
	dm.max_global_id = mpiRanks == 1 ? numPoints : GlobalUniqueID(global_max_coords,dm/*.spacing,dm.global_extent,dm.global_bounds*/);	

	ZERO_ID = dm.max_global_id + 1; // this is the of the singularity vector
	/// NOTICE: we do not have to perturbate the singularity vector
	// if(pertubate) {
	// 	double zero_vec[3] = {0,0,0};
	// 	Perturbate(zero_vec, ZERO_ID, max_global_id);		
	// 	singularity[0] += zero_vec[0];
	// 	singularity[1] += zero_vec[1];
	// 	singularity[2] += zero_vec[2];
	// }

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
		//vecCellIds.reserve(numCells * 2);
		vecCellIds = new vtkIdType[numCells * 6];
		numCellIds=3;
		numSimplices = numCells * 2;
		numSimplicesPerCell = 2;
	}
	else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType) {
		//vecCellIds.reserve(numCells * 5);
		vecCellIds = new vtkIdType[numCells * 20];
		numCellIds=4;
		numSimplices = numCells * 5;
		numSimplicesPerCell = 5;
	}
	else if (VTK_TRIANGLE == cellType) {
		//vecCellIds.reserve(numCells);
		vecCellIds = new vtkIdType[numCells * 3];
		numCellIds=3;
		numSimplices = numCells;
		numSimplicesPerCell = 1;
	}
	else if (VTK_TETRA == cellType) {
		//vecCellIds.reserve(numCells);
		vecCellIds = new vtkIdType[numCells * 4];
		numCellIds=4;
		numSimplices = numCells;
		numSimplicesPerCell = 1;
	}	

	#pragma omp parallel firstprivate(vecCellPerThread) //private(dm,mpiRanks)
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();//Thread ID
		double v[3], p[3]; 

		/// NOTICE: the following two function initialize the two arrays following a spatial-aware strategy.
		///		That said, since we use a vtkDataSet object, this is already guaranteed in the input.
		///		Then, the performance are equivalent, and, thus, we kept the simpler method.
		/*if(VTK_PIXEL == cellType || VTK_QUAD == cellType)
			InitializePointsArray_2D(input,vectors,dm,mpiRanks);
		else if(VTK_VOXEL == cellType || VTK_TETRA == cellType)
			InitializePointsArray_3D(input,vectors,dm,mpiRanks);
		else */
		{
			 #pragma omp for nowait
			 for(vtkIdType i=0; i < numPoints; i++) 
			 {								
			 	vectors->GetTuple(i,v);//&vector[i * 3]);
			 	input->GetPoint(i,p);//&position[i * 3]);
			
				vector[i*3] = static_cast<mpf_float_custom>(v[0]);
				vector[i*3+1] = static_cast<mpf_float_custom>(v[1]);
				vector[i*3+2] = static_cast<mpf_float_custom>(v[2]);

				position[i*3] = static_cast<mpf_float_custom>(p[0]);
				position[i*3+1] = static_cast<mpf_float_custom>(p[1]);
				position[i*3+2] = static_cast<mpf_float_custom>(p[2]);

				 // -- if we have just one MPI process then we can directly use the point id, since the indexing is given and consistent
			 	// -- otherwise, in case of multiple MPI processes we have to derive the global id of the point from some geometric information linked to the grid
			 	vtkIdType global_id = mpiRanks == 1 ? i : GlobalUniqueID(&position[i * 3],dm);						 		
			 	Perturbate(&vector[i * 3], global_id, dm.max_global_id);
			 }
			//std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Error: bad cell type! MPI is not supported here " << std::endl;			
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
				if(i%2){
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
				std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Error: unknown cell type " << std::endl;
				continue;
			}	
			
		}
	}

	max_memory += numCells*numSimplicesPerCell*numCellIds*sizeof(vtkIdType);
	if(mpiRank == 0) std::cout<<"Memory read/write 3D-case: "<<max_memory<<" (bytes) "<<max_memory / std::pow(1024,2) << "(MBs)"<<std::endl;	

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
	//if(mpiRank == 0) 
	std::cout << "[MPI:" << mpiRank << "] [CriticalPointExtractor::identify_critical_points] Extracted " << numSimplices/*vecCellIds.size()*/ << " simplices" << std::endl;
}

/*void CriticalPointExtractor::InitializePointsArray_2D(vtkDataSet * input, vtkDataArray * vectors, DataSetMetadata &dm, int &mpiRanks) 
{
	int j_pos = -1, i_pos = -1; 					
	if (iExchangeIndex == 0) {j_pos =2; i_pos = 1;} 	//2D dataset with yz
	else if (iExchangeIndex == 1) {j_pos = 2; i_pos = 0;} 	//2D dataset with xz
	else if (iExchangeIndex == 2) {j_pos = 1; i_pos = 0;}	//2D dataset with xy
    
#ifdef _WIN32
	#pragma omp for nowait
#else
	#pragma omp for collapse(2)
#endif
	for(vtkIdType j=0; j < dm.dimensions[j_pos]; j+=2) { //Y				
		for(vtkIdType i=0; i < dm.dimensions[i_pos]; i+=2) { //X
				
			/// these 2 loops iterate over one line					
			for(vtkIdType y=0; y < 2; y++) {
				for(vtkIdType x=0; x < 2; x++) {
					vtkIdType id = i + x + (j+y)*dm.dimensions[i_pos];
					vectors->GetTuple(id,&vector[id * 3]);
					input->GetPoint(id, &position[id * 3]);

					// -- if we have just one MPI process then we can directly use the point id, since the indexing is given and consistent
					// -- otherwise, in case of multiple MPI processes we have to derive the global id of the point from some geometric information linked to the grid
					vtkIdType global_id = mpiRanks == 1 ? id : GlobalUniqueID(&position[id * 3],dm);			
					Perturbate(&vector[id * 3], global_id, dm.max_global_id);
					//touched[id]++;					
				}	
			}
		}
	}
}

void CriticalPointExtractor::InitializePointsArray_3D(vtkDataSet * input, vtkDataArray * vectors, DataSetMetadata &dm, int &mpiRanks) 
{	
#ifdef _WIN32
	#pragma omp for nowait
#else
	#pragma omp for collapse(3)
#endif
	for(vtkIdType w=0; w < dm.dimensions[2]; w+=2) { //Z
		for(vtkIdType j=0; j < dm.dimensions[1]; j+=2) { //Y				
			for(vtkIdType i=0; i < dm.dimensions[0]; i+=2) { //X
				
				/// these 3 loops iterate over one line					
				for(vtkIdType z=0; z < 2; z++) {
					for(vtkIdType y=0; y < 2; y++) {
						for(vtkIdType x=0; x < 2; x++) {
							vtkIdType id = i + x + (j+y)*dm.dimensions[0] + (w+z)*dm.dimensions[0]*dm.dimensions[1];
							vectors->GetTuple(id,&vector[id * 3]);
							input->GetPoint(id, &position[id * 3]);
							// -- if we have just one MPI process then we can directly use the point id, since the indexing is given and consistent
							// -- otherwise, in case of multiple MPI processes we have to derive the global id of the point from some geometric information linked to the grid
							vtkIdType global_id = mpiRanks == 1 ? id : GlobalUniqueID(&position[id * 3],dm);			
							Perturbate(&vector[id * 3], global_id, dm.max_global_id);								
						}	
					}			
				}
			}
		}
	}	
}*/

vtkIdType CriticalPointExtractor::GlobalUniqueID(mpf_float_custom* pos, DataSetMetadata &dm/*double *spacing, int *global_extent, double * global_bounds*/)
{
	///Function that calculates global unique id

	/// 1. transpose the point coordinates to the positive range
	double posX = pos[0].convert_to<double>() - dm.global_bounds[0];
	double posY = pos[1].convert_to<double>() - dm.global_bounds[1];
	double posZ = pos[2].convert_to<double>() - dm.global_bounds[2];
	
	/// 2. Compute structured coordinates
	long x = std::lround(posX/dm.spacing[0]);
	long y = std::lround(posY/dm.spacing[1]);
	long z = std::lround(posZ/dm.spacing[2]);

	/// 3. then compute the resolution
	long resx = dm.global_extent[1]+1;
	long resy = dm.global_extent[3]+1;
	long resz = dm.global_extent[5]+1;

	/// 4. then the global id
	// z * xDim * yDim + y * zDim + x
	vtkIdType globalid = z * resx * resy + y * resz + x;

	return globalid;
}

void CriticalPointExtractor::Perturbate(mpf_float_custom* values, vtkIdType &id, vtkIdType &max_global_id) {
	// perturbation function f(e,i,j) = eps^2^i*delta-j
	// eps ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	//eps and delta are constant.. so I compute them one time at the beginning
	vtkIdType i = id + 1;
	double i_norm = static_cast<double>(i)/static_cast<double>(max_global_id);
	double exp_coeff = i_norm*delta;
	// mpf_float_custom exp_coeff = i*delta;
	// double j_norm;

	// if(id==35850 || id==37130 || id==37131) {		
	// 	std::cout<<"eps: "<<eps<<std::endl;
	// 	std::cout<<"Point "<<id<<" values-orig: "<<values[0]<<" "<<values[1]<<" "<<values[2]<<std::endl;
	// }
	// mpf_float_custom exp;// mpf_init(exp);
	// mpf_float_custom const2; const2 = 2; //mpf_init_set_si(const2,2);

	for(int j=0; j<3; j++) {
		// j_norm = static_cast<double>(j+1)/3; //since we are normalizing the point id, we need to normalize as well the j-id --> to keep the perturbation small
		values[j] += boost::multiprecision::pow(eps,std::pow(2,exp_coeff-j_norm[j]));

		// values[j] += boost::multiprecision::pow(eps,boost::multiprecision::pow(2,exp_coeff-(j+1)));

		// mpf_float_custom meps = eps;
		// mpf_float_custom res = boost::multiprecision::pow(meps,std::pow(2,exp_coeff-j_norm));
		// exp = std::pow(2,exp_coeff-j);
		// values[j] = std::pow(eps,exp.get_d());
	}

	// if(id==35850 || id==37130 || id==37131) {		
	// 	std::cout<<"Point "<<id<<" values-pert: "<<values[0]<<" "<<values[1]<<" "<<values[2]<<std::endl;
	// }
	/// FOR DEBUG ONLY --> a perturbation should never be 0
	// if(values[0] == 0 || values[1] == 0 || values[2] == 0) {
	// 	std::cout << "i_norm on id: " << id << " i_norm " << i_norm << " " << exp_coeff << std::endl;
	// 	std::cout << "perturbation on id: " << id << " " << values[0] << " " << values[1] << " " << values[2] << std::endl;
	// }
}

//----------------------------------------------------------------------------
void CriticalPointExtractor::ComputeCriticalCells() 
{
	//vtkIdType cells_num = vecCellIds.size();
	
	int matrixSize = 3;
	if(iExchangeIndex == 3)
		matrixSize = 4;

    std::vector<DynamicMatrix> vecMatrices;
    vecMatrices.resize(numThreads);
	if(matrixSize == 4)
		vecMatrices.assign(numThreads,DynamicMatrix(4,4));//Eigen::Matrix4d());
	if(matrixSize == 3)
		vecMatrices.assign(numThreads,DynamicMatrix(3,3));//Eigen::Matrix3d());	

	//Check for every cell if a critical point (passed singularity as argument) exists
#pragma omp parallel firstprivate(vecMatrices) //private(vector,position)//private(singularity,ZERO_ID)
	{
		//Local variables per thread
		int threadIdx = omp_get_thread_num();	//Thread ID

		//Remove synchronization with nowait
		#pragma omp for
		for (vtkIdType i = 0; i < numSimplices; i++) {
			//If the cell contains a the singularity add them to the output and we can break
			CriticalPointType ret = PointInCell(&vecCellIds[i*numCellIds], vecMatrices[threadIdx]);
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
	singularityType->SetName("singularityType");

	//Prepare output data sequentially. Insert every critical cell to output
	for (auto cellID : vecCriticalCellIDs)
	{
		vtkIdType positionInArray = cellID.id * numCellIds;
		const vtkIdType* vecVertexIds = &vecCellIds[positionInArray];
		vtkSmartPointer<vtkIdList> newPointIDs = vtkSmartPointer<vtkIdList>::New();	
		
		if(numCellIds == 4)
		{
			//Insert the type of singularity 
			singularityType->InsertNextTuple1(cellID.type);

			double pCoords1[3];
			pCoords1[0] = position[vecVertexIds[0] * 3].convert_to<double>();
			pCoords1[1] = position[vecVertexIds[0] * 3 + 1].convert_to<double>();
			pCoords1[2] = position[vecVertexIds[0] * 3 + 2].convert_to<double>();
			
			double pCoords2[3];
			pCoords2[0] = position[vecVertexIds[1] * 3].convert_to<double>();
			pCoords2[1] = position[vecVertexIds[1] * 3 + 1].convert_to<double>();
			pCoords2[2] = position[vecVertexIds[1] * 3 + 2].convert_to<double>();
			
			double pCoords3[3];
			pCoords3[0] = position[vecVertexIds[2] * 3].convert_to<double>();
			pCoords3[1] = position[vecVertexIds[2] * 3 + 1].convert_to<double>();
			pCoords3[2] = position[vecVertexIds[2] * 3 + 2].convert_to<double>();

			double pCoords4[3];
			pCoords4[0] = position[vecVertexIds[3] * 3].convert_to<double>();
			pCoords4[1] = position[vecVertexIds[3] * 3 + 1].convert_to<double>();
			pCoords4[2] = position[vecVertexIds[3] * 3 + 2].convert_to<double>();

			double midPoint[3];
			midPoint[0] = (pCoords1[0] + pCoords2[0] + pCoords3[0] + pCoords4[0]) / 4;
			midPoint[1] = (pCoords1[1] + pCoords2[1] + pCoords3[1] + pCoords4[1]) / 4;
			midPoint[2] = (pCoords1[2] + pCoords2[2] + pCoords3[2] + pCoords4[2]) / 4;

			vtkIdType newPointID = pointArray->InsertNextPoint(midPoint[0], midPoint[1], midPoint[2]);
			newPointIDs->InsertNextId(newPointID);
		}
		if(numCellIds == 3)
		{
			//Insert the type of singularity 
			singularityType->InsertNextTuple1(cellID.type);

			double pCoords1[3];
			pCoords1[0] = position[vecVertexIds[0] * 3].convert_to<double>();
			pCoords1[1] = position[vecVertexIds[0] * 3 + 1].convert_to<double>();
			pCoords1[2] = position[vecVertexIds[0] * 3 + 2].convert_to<double>();
			
			double pCoords2[3];
			pCoords2[0] = position[vecVertexIds[1] * 3].convert_to<double>();
			pCoords2[1] = position[vecVertexIds[1] * 3 + 1].convert_to<double>();
			pCoords2[2] = position[vecVertexIds[1] * 3 + 2].convert_to<double>();
			
			double pCoords3[3];
			pCoords3[0] = position[vecVertexIds[2] * 3].convert_to<double>();
			pCoords3[1] = position[vecVertexIds[2] * 3 + 1].convert_to<double>();
			pCoords3[2] = position[vecVertexIds[2] * 3 + 2].convert_to<double>();

			double midPoint[3];
			midPoint[0] = (pCoords1[0] + pCoords2[0] + pCoords3[0] ) / 3; 
			midPoint[1] = (pCoords1[1] + pCoords2[1] + pCoords3[1] ) / 3; 
			midPoint[2] = (pCoords1[2] + pCoords2[2] + pCoords3[2] ) / 3; 
			
			vtkIdType newPointID = pointArray->InsertNextPoint(midPoint[0], midPoint[1], midPoint[2]);
			newPointIDs->InsertNextId(newPointID);
		}
		cellArray->InsertNextCell(newPointIDs);
	}
		
	//Add points and cells to polydata
	outputData->SetPoints(pointArray); 
	outputData->SetCells(VTK_VERTEX, cellArray);
	outputData->GetPointData()->AddArray(singularityType);
}

CriticalPointExtractor::CriticalPointType CriticalPointExtractor::PointInCell(const vtkIdType* ids, DynamicMatrix &vecMatrix) {		
	std::array<vtkIdType, 4> tmpIds;
	int numIds = numCellIds;
	if(numIds == 3)
		tmpIds = {ids[0],ids[1],ids[2]};
	else
		tmpIds = {ids[0],ids[1],ids[2],ids[3]};

	// 1. compute the initial sign of the determinant of the cell
	mpf_float_custom targetDeterminant = ComputeDeterminant(tmpIds, vecMatrix, false, 0);
	bool targetDirection = DeterminantCounterClockWise(targetDeterminant);
	//Check for non data values (vector is zero and determinant also) 
	if (targetDeterminant == 0)
	{
		this->deg_cases++;
		return REGULAR;
	}
	
	mpf_float_custom tmpDeterminant;
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
			return REGULAR; // regular cell
		}
	}

	if(numIds == 3)
		tmpIds = {ids[0],ids[1],ids[2]};
	else
		tmpIds = {ids[0],ids[1],ids[2],ids[3]};
	mpf_float_custom initialDeterminant = ComputeDeterminant(tmpIds, vecMatrix, true);	
	bool initialDirection     = DeterminantCounterClockWise(initialDeterminant);	

	if (initialDirection != targetDirection)
	{
		return SADDLE; // we found a saddle
	}

	return SINGULARITY; // the cell is critical, since the sign never change
}

mpf_float_custom CriticalPointExtractor::ComputeDeterminant(	
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

	bool vector_field = false;

	// double vecValues[3];
	for (vtkIdType i = 0; i < numIds; i++) 
	{
		vtkIdType& pointID = tmpIds[i];
		if (pointID != ZERO_ID)
		{
			if(!usePoints)
			{
				vecMatrix(i,0) = vector[pointID*3];//vecVectors[pointID][0];// + vecPerturbation[pointID][0] ;
				vecMatrix(i,1) = vector[pointID*3+1];//vecVectors[pointID][1];// + vecPerturbation[pointID][1] ;
				vecMatrix(i,2) = vector[pointID*3+2];//vecVectors[pointID][2];// + vecPerturbation[pointID][2] ;

				vector_field = true;
			}
			else
			{
				vecMatrix(i,0) = position[pointID*3];//vecPointCoordinates[pointID][0];
				vecMatrix(i,1) = position[pointID*3+1];//vecPointCoordinates[pointID][1];
				vecMatrix(i,2) = position[pointID*3+2];//vecPointCoordinates[pointID][2];
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

	int is_simplex = 0;

	// 2. compute determinant sign
	mpf_float_custom det = 0;
	if(numIds == 4) {		
		det = /*static_cast<Eigen::Matrix4d>*/(vecMatrix).determinant();
		is_simplex = 1;
	}
	else {
		det = /*static_cast<Eigen::Matrix3d>*/(vecMatrix).determinant();	
		is_simplex = 2;
	}

	// if(det == 0) { // special case: all vectors are originally 0, and after the perturbation have small numbers causing an underflow
	// 	vecMatrix = vecMatrix*std::pow(10,28);
	// 	det = static_cast<Eigen::Matrix3d>(vecMatrix).determinant();	
	// }	

	// 3. check the number of swap operation while sorting
	if (swapOperations % 2 != 0) //Odd
	{
		det *= -1;
	}

	// if(det == 0) {
	// 	std::cout<<"are we using vector fields? "<<vector_field<<std::endl;
	// 	std::cout<<"which simplex? "<<is_simplex<<std::endl;
	// 	std::cout<<tmpIds[0]<<" "<<tmpIds[1]<<" "<<tmpIds[2]<<" "<<tmpIds[3]<<std::endl;		
	// 	std::cout<<static_cast<double>(tmpIds[0])/40960<<" "<<static_cast<double>(tmpIds[1])/40960<<" "<<static_cast<double>(tmpIds[2])/40960<<" "<<static_cast<double>(tmpIds[3])/40960<<std::endl;
	// 	std::cout<<vecMatrix<<std::endl;
	// 	std::cout<<det<<std::endl;
	// 	vecMatrix = vecMatrix*std::pow(10,14);
	// 	std::cout<<vecMatrix<<std::endl;
	// 	std::cout<<static_cast<Eigen::Matrix4d>(vecMatrix).determinant()<<std::endl;
	// 	// Eigen::SparseLU<Eigen::Matrix4d>   solver(static_cast<Eigen::Matrix4d>(vecMatrix));
	// 	// solver.logAbsDeterminant();
	// 	int a; std::cin>>a;
	// }

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

bool CriticalPointExtractor::DeterminantCounterClockWise(mpf_float_custom& det)
{
	if (det > 0)
		return true;
	else
		return false;
}
