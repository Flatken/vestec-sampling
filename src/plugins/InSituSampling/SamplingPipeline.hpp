
#include "CPTimeSource.hpp"

#include <vtkCPPipeline.h>
#include <vtkDataObject.h>
#include <vtkCPDataDescription.h>
#include <vtkPStreaklineFilter.h>
#include <vtkPolyData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vector>
#include <map>

class SamplingPipeline : public vtkCPPipeline
{
public: 
	static SamplingPipeline *New();
	SamplingPipeline();
	~SamplingPipeline();

  void InitTimingInformation(int numTimeSteps, double deltaT);
	//Override
	int RequestDataDescription(vtkCPDataDescription *dataDescription) override;
	int CoProcess(vtkCPDataDescription *dataDescription) override;
	int Finalize() override;
private:
	CPTimeSource * m_pInput;
	vtkDataObject * m_pCurrentGrid;
	vtkDataObject * m_pLastGrid;
	vtkPStreaklineFilter* m_pStreakLineFilter;
	vtkPolyData* m_pSeeds;
	std::vector<double> timesteps;
  double deltaT;
	double previousTime;
};
