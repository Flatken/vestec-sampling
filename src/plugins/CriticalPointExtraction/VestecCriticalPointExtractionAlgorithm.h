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


//typedef Eigen::Matrix<double, 3, 3> Matrix33;
//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DynamicMatrix;
typedef Eigen::Matrix<double, 3, 3> DynamicMatrix;

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
    CriticalPointExtractor() {}
    void identify_critical_points(	vtkSmartPointer<vtkDataSet> input,
                      vtkSmartPointer<vtkDataSet> output, std::vector<double*> singlarities);

  private:
	int  Sort(vtkSmartPointer<vtkIdList> ids);
    int	 Sort3(vtkSmartPointer<vtkIdList> ids);
	int  Sort4(vtkSmartPointer<vtkIdList> ids);
    double Positive(vtkSmartPointer<vtkIdList> ids, vtkSmartPointer<vtkDataSet> grid, double *currentSingularity, long pertubationID = -1);
    bool PointInCell(vtkSmartPointer<vtkIdList> ids, vtkSmartPointer<vtkDataSet> grid, double* currentSingularity);
    double toFixed(double val);
    bool DeterminatCounterClockWise(double det);
	vtkIdType ZERO_ID;
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