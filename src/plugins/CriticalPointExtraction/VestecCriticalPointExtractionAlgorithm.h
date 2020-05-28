#ifndef __VESTECSEEDINGALGORITHM_H
#define __VESTECSEEDINGALGORITHM_H

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPStreaklineFilter.h>
#include <vtkPParticleTracer.h>
//#include <vtkCell3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkTetra.h>

#include <map>
#include <string>
#include <random>
class vtkDataSet;


// class AbstractGenerator
// {
//   public:
//     virtual ~AbstractGenerator(){};
//     virtual double gen() = 0;
// };

// // A template as generic number generator
// template<class T> class RandomNumberGenerator : public AbstractGenerator
// {
// public:
//     RandomNumberGenerator(double min, double max) 
//     : generator(rd()),
//       distribution(min, max)
//     {}
//     ~RandomNumberGenerator() override {}

//     double gen() override
//     {
//       return distribution(generator);
//     } 
// private:
//     double min;
//     double max;

//     std::random_device rd;
//     std::mt19937 generator; 
//     T distribution;
// };

/// This class implements the algorithm described in:
/// "Robust Detection of Singularities in Vector Fields" by Bhatia et al.
/// -----
/// add other references here:
/// 
class CriticalPointExtractor {
  public:
    CriticalPointExtractor() {}
    vtkSmartPointer<vtkPolyData> identify_critical_points(vtkSmartPointer<vtkDataSet> grid);

  private:
    // float compute_degree(vtkTetra* tet);
    // float compute_solid_angle(vtkTetra* tet);
    int Positive(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid);
    bool PointInCell(vtkCell *cell, vtkSmartPointer<vtkDataSet> grid);
};


class VTK_EXPORT VestecCriticalPointExtractionAlgorithm : public vtkPolyDataAlgorithm
{
public:
  static VestecCriticalPointExtractionAlgorithm *New();
  vtkTypeMacro(VestecCriticalPointExtractionAlgorithm,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkDataSet* GetOutput();
  vtkDataSet* GetOutput(int);
  virtual void SetOutput(vtkDataObject* d);

  // Description:
  // see vtkAlgorithm for details
  virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*) override;

  // this method is not recommended for use, but lots of old style filters use it
  vtkDataObject* GetInput();
  vtkDataObject* GetInput(int port);
  vtkDataSet* GetLabelHierarchyInput(int port);

  // Description:
  // Set an input of this algorithm. You should not override these
  // methods because they are not the only way to connect a pipeline.
  // Note that these methods support old-style pipeline connections.
  // When writing new code you should use the more general
  // vtkAlgorithm::SetInputConnection().  These methods transform the
  // input index to the input port index, not an index of a connection
  // within a single port.
  void SetInput( vtkDataObject* );
  void SetInput( int, vtkDataObject* );

  // Description:
  // Add an input of this algorithm.  Note that these methods support
  // old-style pipeline connections.  When writing new code you should
  // use the more general vtkAlgorithm::AddInputConnection().  See
  // SetInput() for details.
  void AddInput( vtkDataObject* );
  void AddInput( int, vtkDataObject* );

  /**
   * Set the number of created seeds around an input point
   */
  vtkSetMacro(NumberOfPointsAroundSeed, int);
  vtkGetMacro(NumberOfPointsAroundSeed, int);

  /**
   * Set seeding radius in percent of total domain extents
   */
  vtkSetMacro(PercentOfDomain, double);
  vtkGetMacro(PercentOfDomain, double);

  /**
   * Define the mode how ranRandomNumberGeneratordom seeds are generated within a defined range
   */
  vtkSetMacro(DistributionMode, int);
  vtkGetMacro(DistributionMode, int)

  /**
   * Use the scalar field range to determine the number of seeds 
   */
  vtkSetMacro(UseScalarRange, bool);
  vtkGetMacro(UseScalarRange, bool)

protected:
  VestecCriticalPointExtractionAlgorithm();
  ~VestecCriticalPointExtractionAlgorithm();

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestDataObject(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector );

  // convenience method
  virtual int RequestInformation(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector );

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector );

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestUpdateExtent(
    vtkInformation*,
    vtkInformationVector**,
    vtkInformationVector* );

  virtual int FillOutputPortInformation( int port, vtkInformation* info ) override;
  virtual int FillInputPortInformation( int port, vtkInformation* info ) override;

private:
  VestecCriticalPointExtractionAlgorithm( const VestecCriticalPointExtractionAlgorithm& ); // Not implemented.
  void operator = ( const VestecCriticalPointExtractionAlgorithm& );  // Not implemented.

  int       NumberOfPointsAroundSeed = 10;
  double    PercentOfDomain = 5;

  int DistributionMode = 0;

  bool UseScalarRange = 0;
};


#endif