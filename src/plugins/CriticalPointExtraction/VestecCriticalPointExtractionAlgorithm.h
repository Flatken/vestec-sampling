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


//Matrix to compute the determinat
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DynamicMatrix;

//Forward declaration of vtkDataSet
class vtkDataSet;

/// This class implements the algorithm described in:
/// "Robust Detection of Singularities in Vector Fields" by Bhatia et al.
/// -----
/// add other references here:
/// (extending this method) "Detection and classification of critical points in piecewise linear vector fields" by Wang et al.
/// (survey) "A survey of topology‐based methods in visualization" by Heine et al.
/// (survey) "From numerics to combinatorics: a survey of topological methods for vector field visualization" by Wang et al.
class CriticalPointExtractor {
  public:
    /**
     * Store vector and points in internal datat structure 
     */
    CriticalPointExtractor(vtkSmartPointer<vtkDataSet> input, bool pertubate, double* currentSingularity);

    /**
     * Identify the critical cells 
     */
    void ComputeCriticalCells(vtkSmartPointer<vtkDataSet> output);

    void CleanDuplicates(vtkSmartPointer<vtkPolyData> output);

    enum CriticalPointType { REGULAR=0, SADDLE=-1, SINGULARITY=1 };
    struct CriticalPoint {
      vtkIdType id;
      CriticalPointType type;

      CriticalPoint(vtkIdType i, CriticalPointType t) { 
        id = i; 
        type = t;
      }
    };
    
  private:
    /**
     * Sort the vector integers but returns the needed swap operations
     */
    int  Sort(std::vector<vtkIdType> &ids);

    /**
     * Sort the vector with 3 integers but returns the needed swap operations
     */
    int	 Sort3(std::vector<vtkIdType> &ids);

    /**
     * Sort the vector with 4 integers but returns the needed swap operations
     */
	  int  Sort4(std::vector<vtkIdType> &ids);
    
    /**
     * Compute the determinant
     * tmpIds: Contains the vertex indices needed to fill the matrix. Is a copy since they will be sorted
     * grid: Gives access to vector field
     * currentSingularity: The singularity e.g. zero vector (0,0,0)
     * vecMatrix: The matrix used to compute the determinant
     */
    double ComputeDeterminant(std::vector<vtkIdType> tmpIds, DynamicMatrix &vecMatrix, bool usePoints, long pertubationID = -1);
   
    /**
     * Check if the sigularity is in cell
     * ids: The vertex ids spanning the cell
     * vecMatrix: The matrix used to compute the determinant
     */ 
    CriticalPointType PointInCell(std::vector<vtkIdType> &ids, DynamicMatrix &vecMatrix);

    /**
     * Check if direction of the determinat is positive (counter-clockwise) 
     */
    bool DeterminatCounterClockWise(double& det);

    /**
     * Double to fixed precision and pertubation based on id
     */
    void Pertubate(double* values, vtkIdType id);

    
private:
    vtkIdType ZERO_ID;  //!< Vertex ID of the singularity: always number of vertices + 1 
    int iExchangeIndex; //!< The row id of the matrix to exchange with the singularity    
    std::vector<double*> vecPointCoordinates; //!< Store point coordinates
    std::vector<double*> vecVectors; //!< Store vector field
    std::vector<std::vector<vtkIdType>> vecCellIds;  //!< The point ids for each cell
    double singularity[3]; //!< The singularity to identify
    int numThreads = 12; //!< Number of OpenMP threads
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