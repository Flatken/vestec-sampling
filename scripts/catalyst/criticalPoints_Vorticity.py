
import sys
print("The Python version is %s.%s.%s" % sys.version_info[:3])

from paraview.simple import *
from paraview import coprocessing


# --------------------------------------------------------------
# The following loads TTK's plugins.
# Topology Toolkit 0.9.7
# --------------------------------------------------------------
import glob
import os
from os.path import join as ttk_path_join

ttk_plugins_path = "../lib/plugins/"
for x in glob.glob(
    ttk_path_join(ttk_plugins_path, "*.so" if os.name == "posix" else "*.dll")
):
    LoadPlugin(x, ns=globals())
    print("Loading Plugin: "+x)
#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.5.2

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.5.2

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.5.2

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      grid_ = coprocessor.CreateProducer(datadescription, 'input')
      
      # create a new 'Compute Derivatives'
      computeVorticity = ComputeDerivatives(Input=grid_)
      computeVorticity.Scalars = ['POINTS', 'p']
      computeVorticity.Vectors = ['POINTS', 'Vec']
      computeVorticity.OutputVectorType = 'Vorticity'
      computeVorticity.OutputTensorType = 'Nothing'

      # create a new 'Calculator'
      vorticityMagnitutde = Calculator(Input=computeVorticity)
      vorticityMagnitutde.AttributeType = 'Cell Data'
      vorticityMagnitutde.ResultArrayName = 'vortMag'
      vorticityMagnitutde.Function = 'Vorticity_Z'

      # create a new 'Cell Data to Point Data'
      grid = CellDatatoPointData(Input=vorticityMagnitutde)

      # create a new 'Tetrahedralize'
      tetrahedralize1 = Tetrahedralize(Input=grid)
      
      # create a new 'TTK PersistenceCurve'
      tTKPersistenceCurve1 = TTKPersistenceCurve(Input=tetrahedralize1)
      tTKPersistenceCurve1.ScalarField = 'vortMag'
      tTKPersistenceCurve1.InputOffsetField = 'vortMag'

      # create a new 'TTK PersistenceDiagram'
      tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=tetrahedralize1)
      tTKPersistenceDiagram1.ScalarField = 'vortMag'
      tTKPersistenceDiagram1.InputOffsetField = 'vortMag'

      # create a new 'Threshold'
      positives = Threshold(Input=tTKPersistenceDiagram1)
      positives.Scalars = ['CELLS', 'PairIdentifier']
      positives.ThresholdRange = [0.01, 9184.0]

      # create a new 'Threshold'
      persistentThreshold = Threshold(Input=positives)
      persistentThreshold.Scalars = ['CELLS', 'Persistence']
      persistentThreshold.ThresholdRange = [1.1, 17.2276411574604]

      # create a new 'Extract Surface'
      extractSurface1 = ExtractSurface(Input=persistentThreshold)

      # create a new 'TTK TopologicalSimplification'
      tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=tetrahedralize1,
          Constraints=persistentThreshold)
      tTKTopologicalSimplification1.ScalarField = 'vortMag'
      tTKTopologicalSimplification1.InputOffsetField = 'vortMag'
      tTKTopologicalSimplification1.OutputOffsetScalarField = ''

      # create a new 'TTK ScalarFieldCriticalPoints'
      tTKScalarFieldCriticalPoints1 = TTKScalarFieldCriticalPoints(Input=tTKTopologicalSimplification1)
      tTKScalarFieldCriticalPoints1.ScalarField = 'vortMag'
      tTKScalarFieldCriticalPoints1.ForceInputOffsetField = 1
      tTKScalarFieldCriticalPoints1.InputOffsetfield = 'ttkOffsetScalarField'
      tTKScalarFieldCriticalPoints1.Withboundarymask = 0

      # create a new 'Threshold'
      diagonal = Threshold(Input=tTKPersistenceDiagram1)
      diagonal.Scalars = ['CELLS', 'PairIdentifier']
      diagonal.ThresholdRange = [-1.0, -0.1]

      # Compute streakline snippets
      vestecSamplingAlgorithm = VestecSamplingAlgorithm(Grid=grid_,Seeds=tTKScalarFieldCriticalPoints1)
      vestecSamplingAlgorithm.Vectors = ['POINTS', 'Vel']
      
      #SetActiveSource(vestecSamplingAlgorithm)

      # Now any catalyst writers
      xMLPPolyDataWriter1 = servermanager.writers.XMLPPolyDataWriter(Input=vestecSamplingAlgorithm)
      coprocessor.RegisterWriter(xMLPPolyDataWriter1, filename='VestecSamplingAlgorithm1_%t.pvtp', freq=1, paddingamount=0)

      #ParallelUnstructuredGridWriter1 = coprocessor.CreateWriter( XMLPUnstructuredGridWriter, "seeds_%t.pvtu", 1 )

      # ----------------------------------------------------------------
      # finally, restore active source
      # SetActiveSource(grid)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
