#ifndef UnsteadySource_H
#define UnsteadySource_H

#include <vtkAlgorithm.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <map>
#include <vector>

using namespace std;

/**
 * Helper class to define a time dependent data source for our catalyst pipeline
 */
class UnsteadySource : public vtkAlgorithm
{
public:
  static UnsteadySource *New();
  vtkTypeMacro(UnsteadySource,vtkAlgorithm);

  void SetBoundingBox(double x0, double x1, double y0,
                      double y1, double z0, double z1)
  {
    this->BoundingBox[0] = x0;
    this->BoundingBox[1] = x1;
    this->BoundingBox[2] = y0;
    this->BoundingBox[3] = y1;
    this->BoundingBox[4] = z0;
    this->BoundingBox[5] = z1;
  }

  void SetExtent(int xMin, int xMax, int yMin, int yMax,
                    int zMin, int zMax)
  {
    int modified = 0;

    if (this->Extent[0] != xMin)
    {
      modified = 1;
      this->Extent[0] = xMin ;
    }
    if (this->Extent[1] != xMax)
    {
      modified = 1;
      this->Extent[1] = xMax ;
    }
    if (this->Extent[2] != yMin)
    {
      modified = 1;
      this->Extent[2] = yMin ;
    }
    if (this->Extent[3] != yMax)
    {
      modified = 1;
      this->Extent[3] = yMax ;
    }
    if (this->Extent[4] != zMin)
    {
      modified = 1;
      this->Extent[4] = zMin ;
    }
    if (this->Extent[5] != zMax)
    {
      modified = 1;
      this->Extent[5] = zMax ;
    }
    if (modified)
    {
      this->Modified();
    }
  }

   int GetNumberOfTimeSteps()
   {
     return static_cast<int>(this->TimeSteps.size());
   }

   void SetNumberOfTimeSteps(std::vector<double>& t)
   {
     this->TimeSteps = t;
     this->GetInformation()->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),&t[0],static_cast<int>(t.size()));
     double range[2]= {0,t.size()};
     this->GetInformation()->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),range,2);
     std::cout << "Timesteps: " << this->TimeSteps.size() << std::endl;
     for (auto &value : this->TimeSteps) // access by reference to avoid copying
       std::cout << " " << value;
     std::cout << std::endl;

    
   }

	void InsertDataObject(double t, vtkDataObject* obj)
   {
		  m_pDataObjects.insert(std::make_pair(t, obj));
   }

   vtkDataObject* GetDataObject(double t)
   {
	   auto result = m_pDataObjects.find(t);
	   if(result != m_pDataObjects.end())
	   		return result->second;
		else
		{
			return nullptr;
		}
		
   }
 
 protected:
   UnsteadySource()
   {
     this->SetNumberOfInputPorts(0);
     this->SetNumberOfOutputPorts(1);
 
     this->Extent[0] = 0;
     this->Extent[1] = 1;
     this->Extent[2] = 0;
     this->Extent[3] = 1;
     this->Extent[4] = 0;
     this->Extent[5] = 1;
 
     this->BoundingBox[0]=0;
     this->BoundingBox[1]=1;
     this->BoundingBox[2]=0;
     this->BoundingBox[3]=1;
     this->BoundingBox[4]=0;
     this->BoundingBox[5]=1;
   }
   ~UnsteadySource() { }

   void GetSpacing(double dx[3])
   {
     for(int i=0; i<3; i++)
     {
       dx[i] = (this->BoundingBox[2*i+1]- this->BoundingBox[2*i]) / (this->Extent[2*i+1] - this->Extent[2*i]);
     }
   }

   int ProcessRequest(vtkInformation* request,
                      vtkInformationVector** inputVector,
                      vtkInformationVector* outputVector) override
   {
     // generate the data
     if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
     {
       return this->RequestData(request, inputVector, outputVector);
     }
 
     // execute information
     if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
     {
       return this->RequestInformation(request, inputVector, outputVector);
     }
     return this->Superclass::ProcessRequest(request, inputVector, outputVector);
   }
 
   int FillOutputPortInformation(int, vtkInformation *info) override
   {
     info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
     return 1;
   }
 
   virtual int RequestInformation(vtkInformation *,
                                  vtkInformationVector **,
                                  vtkInformationVector *outputInfoVector)
   {
     // get the info objects
     vtkInformation *outInfo = outputInfoVector->GetInformationObject(0);
 
     double range[2]= {0,this->TimeSteps.size()};
     outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),range,2);
     outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),&this->TimeSteps[0], static_cast<int>(this->TimeSteps.size()));
     outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), this->Extent,6);
  
     double spacing[3];
     this->GetSpacing(spacing);
 
     outInfo->Set(vtkDataObject::SPACING(), spacing[0], spacing[1], spacing[2]);
 
     double origin[3] = {this->BoundingBox[0],this->BoundingBox[2],this->BoundingBox[4]};
     outInfo->Set(vtkDataObject::ORIGIN(),origin,3);
     outInfo->Set(vtkAlgorithm::CAN_PRODUCE_SUB_EXTENT(), 1);
     return 1;
   }
 
 
   int RequestData(
     vtkInformation* ,
     vtkInformationVector** vtkNotUsed( inputVector ),
     vtkInformationVector* outputVector)
   {
     vtkInformation *outInfo = outputVector->GetInformationObject(0);
     vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());
 
     double timeStep = outInfo->Get( vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
     output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),timeStep);
     // set the extent to be the update extent
     

    std::cout << "Requesting data for timestep: " << timeStep <<  std::endl;
    //TODO Pass correct dataset
    vtkDataObject* internal = GetDataObject(timeStep);
    if(internal)
    {
      output->ShallowCopy(internal);
      return 1;
    }
    else
    {
      std::cout << "Error: Requesting data for timestep: " << timeStep << " not found" << std::endl;
      return 0;
    }
   }

  
 private:
   UnsteadySource(const UnsteadySource&) = delete;
   void operator=(const UnsteadySource&) = delete;
 
   std::vector<double> TimeSteps;
   std::map<double,vtkDataObject*> m_pDataObjects;
   int Extent[6];
   double BoundingBox[6];
   int Spacing;
 
 };

 #endif // UnsteadySource