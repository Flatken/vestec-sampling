#ifndef __VESTECSEEDINGALGORITHM_H
#define __VESTECSEEDINGALGORITHM_H

#include <vtkAlgorithm.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkPolyData.h>
#include <vtkCell.h>

#include <map>
#include <string>
#include <vector>
#include <random>

#include <Eigen/Dense>
#include <cnl/all.h>

//Matrix to compute the determinant
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 1, 4, 4> DynamicMatrix;
//Forward declaration of vtkDataSet
class vtkDataSet;

/// This class extends the algorithm described in:
/// "Detection and classification of critical points in piecewise linear vector fields" by Wang et al.
/// (that extends the following method) "Robust Detection of Singularities in Vector Fields" by Bhatia et al./// 
/// -----
/// other references here:
/// (possible extension using a different index for classification) "Morse set classification and hierarchical refinement using conley index" by Chen G et al.
/// (survey) "A survey of topology‐based methods in visualization" by Heine et al.
/// (survey) "From numerics to combinatorics: a survey of topological methods for vector field visualization" by Wang et al.
class CriticalPointExtractor {
  public:
    /**
     * Store vector and points in internal data structure 
     */
    CriticalPointExtractor(vtkSmartPointer<vtkDataSet> input, double* currentSingularity, int mpiRank, bool pertubate = true);

    /**
     * Identify the critical cells 
     */
    void ComputeCriticalCells();

    void writeCriticalCells(vtkSmartPointer<vtkDataSet> output);

    enum CriticalPointType { REGULAR=0, SADDLE=-1, SINGULARITY=1 };
    struct CriticalPoint {
      vtkIdType id;
      CriticalPointType type;

      CriticalPoint(vtkIdType i, CriticalPointType t) { 
        id = i; 
        type = t;
      }
    };

    ~CriticalPointExtractor() {     
      delete position;
	    delete vector;
	    delete perturbation;
     
     #pragma omp parallel for
     for(vtkIdType x=0; x < vecCellIds.size(); ++x)
        delete vecCellIds[x];

     vecCellIds.clear();
     vecPointCoordinates.clear();
     vecVectors.clear();
     vecPerturbation.clear();
    }
    
  private:
    /**
     * Sort the vector integers but returns the needed swap operations
     */
    int  Sort(vtkIdType* ids, int n);

    /**
     * Sort the vector with 3 integers but returns the needed swap operations
     */
    int	 Sort3(vtkIdType* ids);

    /**
     * Sort the vector with 4 integers but returns the needed swap operations
     */
	  int  Sort4(vtkIdType* ids);
    
    /**
     * Compute the determinant
     * tmpIds: Contains the vertex indices needed to fill the matrix. Is a copy since they will be sorted
     * grid: Gives access to vector field
     * currentSingularity: The singularity e.g. zero vector (0,0,0)
     * vecMatrix: The matrix used to compute the determinant
     */
    double ComputeDeterminant(const vtkIdType* ids, DynamicMatrix &vecMatrix, bool usePoints, long perturbationID = -1);
   
    /**
     * Check if the singularity is in cell
     * ids: The vertex ids spanning the cell
     * vecMatrix: The matrix used to compute the determinant
     */ 
    CriticalPointType PointInCell(const vtkIdType* ids, DynamicMatrix &vecMatrix);

    /**
     * Check if direction of the determinant is positive (counter-clockwise) 
     */
    bool DeterminantCounterClockWise(double& det);

    /**
     * Double to fixed-precision conversion (based on CNL-library)
     */
    void toFixed(double *values);

    /**
     * Perturbation based on point id
     */
    void Perturbate(double* values, long id, long max_global_id);

    /**
     * Calculate global unique id 
     */
    long GlobalUniqueID(double* pos, double *spacing, int *global_extent, double * global_bounds);
    
private:
    vtkIdType ZERO_ID;  //!< Vertex ID of the singularity: always number of vertices + 1 
    int iExchangeIndex; //!< The row id of the matrix to exchange with the singularity    
    std::vector<double*> vecPointCoordinates; //!< Store point coordinates
    std::vector<double*> vecVectors; //!< Store vector field
    std::vector<double*> vecPerturbation; //!< Store vector field perturbation
    double* position;
	  double* vector;
	  double* perturbation;
    std::vector<vtkIdType*> vecCellIds;  //!< The point ids for each cell
    int numCellIds;
    double singularity[3]; //!< The singularity to identify
    int numThreads; //!< Number of OpenMP threads
    double eps = 1 / std::pow(10,14);
	  double delta = 4; // >=n    
	  std::vector<CriticalPoint> vecCriticalCellIDs; //!< Vector of critical cell ids
};


class VTK_EXPORT VestecCriticalPointExtractionAlgorithm : public vtkDataSetAlgorithm
{
public:
  static VestecCriticalPointExtractionAlgorithm *New();
  vtkTypeMacro(VestecCriticalPointExtractionAlgorithm, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  VestecCriticalPointExtractionAlgorithm();
  ~VestecCriticalPointExtractionAlgorithm();

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector );

  virtual int FillOutputPortInformation(int port, vtkInformation* info) override;

private:
  VestecCriticalPointExtractionAlgorithm( const VestecCriticalPointExtractionAlgorithm& ); // Not implemented.
  void operator = ( const VestecCriticalPointExtractionAlgorithm& );  // Not implemented.
};


#endif