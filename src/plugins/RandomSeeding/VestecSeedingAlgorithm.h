#ifndef __VESTECSEEDINGALGORITHM_H
#define __VESTECSEEDINGALGORITHM_H

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPStreaklineFilter.h>
#include <vtkPParticleTracer.h>
#include <map>

class vtkDataSet;

class VTK_EXPORT VestecSeedingAlgorithm : public vtkPolyDataAlgorithm
{
public:
  static VestecSeedingAlgorithm *New();
  vtkTypeMacro(VestecSeedingAlgorithm,vtkPolyDataAlgorithm);
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


protected:
  VestecSeedingAlgorithm();
  ~VestecSeedingAlgorithm();

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
  VestecSeedingAlgorithm( const VestecSeedingAlgorithm& ); // Not implemented.
  void operator = ( const VestecSeedingAlgorithm& );  // Not implemented.

  int       NumberOfPointsAroundSeed = 10;
  double    PercentOfDomain = 5;
};


#endif