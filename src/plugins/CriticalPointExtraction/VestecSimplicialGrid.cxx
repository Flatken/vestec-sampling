#include "VestecSimplicialGrid.h"

SimplicialGrid::SimplicialGrid(vtkIdType numPoints, vtkIdType numCells, vtkIdType cellType) {
    /// RICCARDO: changed from resize to reserve, since this is going to be the maximum memory needed.. but not actually used

    //Allocate memory for encoding points coordinates and vectors
	// vecPointCoordinates.resize(numPoints);
	// vecVectors.resize(numPoints);
    

    //Allocate size for cells which depends on input cell type
	if(VTK_PIXEL == cellType || VTK_QUAD == cellType) {
		// vecCellIds.resize(numCells * 2);
		numCellIds=3;
		cellIds = new vtkIdType[numCells*2*numCellIds];
        numSimplices = numCells*2;
	}
	else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType) {
		// vecCellIds.resize(numCells * 5);        
		numCellIds=4;
        cellIds = new vtkIdType[numCells*5*numCellIds];
        numSimplices = numCells*5;
	}
	else if (VTK_TRIANGLE == cellType) {
		// vecCellIds.resize(numCells);
		numCellIds=3;
		cellIds = new vtkIdType[numCells*numCellIds];
        numSimplices = numCells;		
	}
	else if (VTK_TETRA == cellType) {
		// vecCellIds.resize(numCells);
		numCellIds=4;
		cellIds = new vtkIdType[numCells*numCellIds];
        numSimplices = numCells;
	}
}

void SimplicialGrid::AddSimplex(
    vtkSmartPointer<vtkDataSet> input, vtkIdType &i, vtkIdList* ids, vtkIdType &cellType, vtkIdType &chunk_size, 
    vtkIdType &mpiRanks, double* spacing, int* global_extent, double* global_bounds, vtkIdType &max_global_id) {

    vtkSmartPointer<vtkDataArray> vectors = input->GetPointData()->GetVectors();	
    vtkIdType numPoints;

	vtkIdType local_id = i % chunk_size;

    if (VTK_PIXEL == cellType || VTK_QUAD == cellType)
	{
		// vtkIdType* test = new vtkIdType[6]{ ids->GetId(0) , ids->GetId(1), ids->GetId(2), 
		// 	ids->GetId(1), ids->GetId(3) , ids->GetId(2)};
		// vecCellIds[i * 2] = &test[0];
		// vecCellIds[i * 2 + 1] = &test[3];
		cellIds[local_id*6] = ids->GetId(0); cellIds[local_id*6+1] = ids->GetId(1); cellIds[local_id*6+2] = ids->GetId(2);
		cellIds[local_id*6+3] = ids->GetId(1); cellIds[local_id*6+4] = ids->GetId(3); cellIds[local_id*6+5] = ids->GetId(2);
        numPoints = 4;
	}
	else if (VTK_VOXEL == cellType || VTK_HEXAHEDRON == cellType)
	{
		// vtkIdType* test = new vtkIdType[20]{ ids->GetId(0) , ids->GetId(6), ids->GetId(4), ids->GetId(5), 
		// 	ids->GetId(3) , ids->GetId(5), ids->GetId(7), ids->GetId(6),
		// 	ids->GetId(3) , ids->GetId(1), ids->GetId(5), ids->GetId(0),
		// 	ids->GetId(0) , ids->GetId(3), ids->GetId(2), ids->GetId(6),
		// 	ids->GetId(0) , ids->GetId(6), ids->GetId(3), ids->GetId(5)};
		// vecCellIds[i * 5] = &test[0];
		// vecCellIds[i * 5 + 1] = &test[4];
		// vecCellIds[i * 5 + 2] = &test[8];
		// vecCellIds[i * 5 + 3] = &test[12];
		// vecCellIds[i * 5 + 4] = &test[16];
        
        cellIds[local_id*20] = ids->GetId(0); cellIds[local_id*20+1] = ids->GetId(6); cellIds[local_id*20+2] = ids->GetId(4); cellIds[local_id*20+3] = ids->GetId(5);
        cellIds[local_id*20+4] = ids->GetId(3); cellIds[local_id*20+5] = ids->GetId(5); cellIds[local_id*20+6] = ids->GetId(7); cellIds[local_id*20+7] = ids->GetId(6);
        cellIds[local_id*20+8] = ids->GetId(3); cellIds[local_id*20+9] = ids->GetId(1); cellIds[local_id*20+10] = ids->GetId(5); cellIds[local_id*20+11] = ids->GetId(0);
        cellIds[local_id*20+12] = ids->GetId(0); cellIds[local_id*20+13] = ids->GetId(3); cellIds[local_id*20+14] = ids->GetId(2); cellIds[local_id*20+15] = ids->GetId(6);
        cellIds[local_id*20+16] = ids->GetId(0); cellIds[local_id*20+17] = ids->GetId(6); cellIds[local_id*20+18] = ids->GetId(3); cellIds[local_id*20+19] = ids->GetId(5);

        numPoints = 8;
       
	}
	else if (VTK_TRIANGLE == cellType)
	{
		// vecCellIds[i] = new vtkIdType[3]{ ids->GetId(0) , ids->GetId(1), ids->GetId(2)};
		cellIds[local_id*3] = ids->GetId(0); cellIds[local_id*3+1] = ids->GetId(1); cellIds[local_id*3+2] = ids->GetId(2);
        numPoints = 3;
	}
	else if (VTK_TETRA == cellType)
	{
		// vecCellIds[i] = new vtkIdType[4]{ ids->GetId(0) , ids->GetId(1), ids->GetId(2), ids->GetId(3)};
		cellIds[local_id*4] = ids->GetId(0); cellIds[local_id*4+1] = ids->GetId(1); cellIds[local_id*4+2] = ids->GetId(2); cellIds[local_id*4+3] = ids->GetId(3);
        numPoints = 4;
	}
	else {
		std::cout << "[MPI:] [SimplicialGrid::AddSimplex] Error: unknown cell type " << std::endl;
		return;
        // continue;
	}

    for(vtkIdType j=0; j < numPoints; j++) 
    {
        vtkIdType v_id = ids->GetId(j);        
        if(vecPoints.find(v_id)==vecPoints.end()) {
            double* c = new double[6];
            // double* v = new double[3];
            
		    input->GetPoint(v_id, &c[0]);	
            vectors->GetTuple(v_id,&c[3]);

		    // -- if we have just one MPI process then we can directly use the point id, since the indexing is given and consistent
		    // -- otherwise, in case of multiple MPI processes we have to derive the global id of the point from some geometric information linked to the grid
		    long global_id = mpiRanks == 1 ? v_id : SimplicialGrid::GlobalUniqueID(&c[0],spacing,global_extent,global_bounds);
			
		    // if(pertubate) {								
		    SimplicialGrid::Perturbate(&c[3], global_id, max_global_id);	
		    // }
			
		    vecPoints[v_id] = c;
        }
    }
	// {
	// 	vectors->GetTuple(i,&vector[i * 3]);
	// 	input->GetPoint(i, &position[i * 3]);	

	// 	// -- if we have just one MPI process then we can directly use the point id, since the indexing is given and consistent
	// 	// -- otherwise, in case of multiple MPI processes we have to derive the global id of the point from some geometric information linked to the grid
	// 	long global_id = mpiRanks == 1 ? i : GlobalUniqueID(&position[i * 3],spacing,global_extent,global_bounds);
			
	// 	// if(pertubate) {								
	// 		Perturbate(&vector[i * 3], global_id, max_global_id);	
	// 	// }
			
	// 	vecVectors[i] 		= &vector[i * 3];
	// 	// vecPerturbation[i] 	= &perturbation[i * 3];
	// 	vecPointCoordinates[i] 	= &position[i * 3];
	// }
}

long SimplicialGrid::GlobalUniqueID(double* pos, double *spacing, int *global_extent, double * global_bounds)
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
	// z * xDim * yDim + y * zDim + x
	long globalid = z * resx * resy + y * resz + x;

	return globalid;
}

void SimplicialGrid::Perturbate(double* values, long id, long max_global_id) {
	// perturbation function f(e,i,j) = eps^2^i*delta-j
	// eps ?? --> constant?
	// i = id (in their implementation is id+1)
	// j = to the component of values --> 0,1,2 (in their implementation is the component +1)

	//eps and delta are constant.. so I compute them one time at the beginning
	vtkIdType i = id + 1;
	double i_norm = static_cast<double>(i)/static_cast<double>(max_global_id);
	double exp_coeff = i_norm*DELTA;
	double j_norm;

	for(int j=0; j<3; j++) {
		j_norm = static_cast<double>(j+1)/3; //since we are normalizing the point id, we need to normalize as well the j-id --> to keep the perturbation small
		values[j] += std::pow(EPS,std::pow(2,exp_coeff-j_norm));
	}

	/// FOR DEBUG ONLY --> a perturbation should never be 0
	// if(values[0] == 0 || values[1] == 0 || values[2] == 0) {
	// 	std::cout << "i_norm on id: " << id << " i_norm " << i_norm << " " << exp_coeff << std::endl;
	// 	std::cout << "perturbation on id: " << id << " " << values[0] << " " << values[1] << " " << values[2] << std::endl;
	// }
}
