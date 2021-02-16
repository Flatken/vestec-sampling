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

#include <map>

const double EPS = 1 / std::pow(10,14);
const double DELTA = 4; // >=n   

class SimplicialGrid {
public:
    SimplicialGrid(vtkIdType numPoints, vtkIdType numCells, vtkIdType cellType);

    ~SimplicialGrid() {     
        // delete position;
	    // delete vector;	  
     
    // #pragma omp parallel for
        // for(vtkIdType x=0; x < vecCellIds.size(); ++x)
        int count = 0;
        for(auto it=vecCellIds.begin(); it!=vecCellIds.end(); ++it) {
            if(count%(this->numCellIds+1)==0)
              delete it->second;
            count++;
        }            
        vecCellIds.clear();

        for(auto it=vecPoints.begin(); it!=vecPoints.end(); ++it) {
            delete it->second.first;
            delete it->second.second;
        }            
        vecPoints.clear();
        // vecVectors.clear();
    }    

    void AddSimplex(vtkSmartPointer<vtkDataSet> input, vtkIdType &i, vtkIdList* ids, vtkIdType &cellType, vtkIdType &chunk_size, 
        vtkIdType &mpiRanks, double* spacing, int* global_extent, double* global_bounds, vtkIdType &max_global_id);

    // inline vtkIdType GetSimplicesNum() { return vecCellIds.size(); }
    // inline vtkIdType* GetSimplex(vtkIdType i) { return vecCellIds[i]; }
    inline std::map<int,vtkIdType*>::iterator GetCellIds_Begin() { return vecCellIds.begin(); }
    inline std::map<int,vtkIdType*>::iterator GetCellIds_End() { return vecCellIds.end(); }
    inline int GetNumCellIds() { return numCellIds; }
    inline vtkIdType GetNumCells() { return numSimplices/*vecCellIds.size()*/; }
    inline vtkIdType* GetSimplex(vtkIdType i) { return &cellIds[i*numCellIds]; }
    inline std::pair<double*,double*>& GetPoint(vtkIdType pId) { /*std::cout<<"numPoints: "<<vecPoints.size()<<std::endl;*/ return vecPoints[pId]; }

    /**
     * Perturbation based on point id
     */
    static void Perturbate(double* values, long id, long max_global_id);

    /**
     * Calculate global unique id 
     */
    static long GlobalUniqueID(double* pos, double *spacing, int *global_extent, double * global_bounds);

private:
    /// the first array encodes the coordinates, while the second the vector associated to a point
    std::map<int,std::pair<double*,double*> > vecPoints; //!< Store point coordinates
    // std::map<int,double*> vecVectors; //!< Store vector field
    // std::vector<double*> vecPerturbation; //!< Store vector field perturbation
    // double* position;
	// double* vector;
	// double* perturbation;
    std::map<int,vtkIdType*> vecCellIds;  //!< The point ids for each cell
    vtkIdType* cellIds;
    int numCellIds;    
    vtkIdType numSimplices;
};