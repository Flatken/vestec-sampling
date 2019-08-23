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
#-----------------------------------
# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.6.0

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.6.0
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      politano_B_ = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Compute Derivatives'
      computeDerivatives1 = ComputeDerivatives(Input=politano_B_)
      computeDerivatives1.Scalars = [None, '']
      computeDerivatives1.Vectors = ['POINTS', 'B']
      computeDerivatives1.OutputVectorType = 'Vorticity'
      computeDerivatives1.OutputTensorType = 'Nothing'

      # create a new 'Calculator'
      calculator1 = Calculator(Input=computeDerivatives1)
      calculator1.AttributeType = 'Cell Data'
      calculator1.ResultArrayName = 'mag'
      calculator1.Function = 'mag(Vorticity)'

      # create a new 'Cell Data to Point Data'
      cellDatatoPointData1 = CellDatatoPointData(Input=calculator1)

      # create a new 'Tetrahedralize'
      tetrahedralize1 = Tetrahedralize(Input=cellDatatoPointData1)

      # create a new 'TTK PersistenceDiagram'
      tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=tetrahedralize1)
      tTKPersistenceDiagram1.ScalarField = 'mag'
      tTKPersistenceDiagram1.InputOffsetField = 'mag'

      # create a new 'Threshold'
      threshold1 = Threshold(Input=tTKPersistenceDiagram1)
      threshold1.Scalars = ['CELLS', 'PairIdentifier']
      threshold1.ThresholdRange = [0.1, 28223.0]

      # create a new 'Threshold'
      persistence = Threshold(Input=threshold1)
      persistence.Scalars = ['CELLS', 'Persistence']
      persistence.ThresholdRange = [0.0025, 0.03259960925605887]

      # create a new 'TTK TopologicalSimplification'
      tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=tetrahedralize1,
          Constraints=persistence)
      tTKTopologicalSimplification1.ScalarField = 'mag'
      tTKTopologicalSimplification1.InputOffsetField = 'mag'
      tTKTopologicalSimplification1.Vertexidentifierfield = 'CriticalType'
      tTKTopologicalSimplification1.OutputOffsetScalarField = ''

      # create a new 'TTK ScalarFieldCriticalPoints'
      tTKScalarFieldCriticalPoints1 = TTKScalarFieldCriticalPoints(Input=tTKTopologicalSimplification1)
      tTKScalarFieldCriticalPoints1.ScalarField = 'mag'
      tTKScalarFieldCriticalPoints1.InputOffsetfield = 'mag'

      # create a new 'Mask Points'
      maskPoints1 = MaskPoints(Input=tTKScalarFieldCriticalPoints1)
      maskPoints1.OnRatio = 1
      maskPoints1.MaximumNumberofPoints = 2000000
      maskPoints1.GenerateVertices = 1
      maskPoints1.SingleVertexPerCell = 1

      # create a new 'Threshold'
      threshold2 = Threshold(Input=maskPoints1)
      threshold2.Scalars = ['POINTS', 'CriticalType']

      # create a new 'Threshold'
      threshold3 = Threshold(Input=maskPoints1)
      threshold3.Scalars = ['POINTS', 'CriticalType']
      threshold3.ThresholdRange = [3.0, 3.0]

      # create a new 'Append Datasets'
      appendDatasets1 = AppendDatasets(Input=[threshold3, threshold2])

      # create a new 'TTK SphereFromPoint'
      vestecSeeding = VestecSeedingAlgorithm(Input=appendDatasets1)
      vestecSeeding.SeedingRadius = 5.0
      vestecSeeding.NumberOfSeeds = 4
      vestecSeeding.RandomDistributionMode = 'Normal distribution'
      vestecSeeding.Array = ['POINTS', 'p']

      # Compute streakline snippets
      vestecSamplingAlgorithm = VestecSamplingAlgorithm(Grid=politano_B_,Seeds=vestecSeeding)
      vestecSamplingAlgorithm.Vectors = ['POINTS', 'B']
      vestecSamplingAlgorithm.IntegrationDuration = 6
      vestecSamplingAlgorithm.StepSize = 1

      SetActiveSource(vestecSamplingAlgorithm)

      # Now any catalyst writers
      xMLPPolyDataWriter1 = servermanager.writers.XMLPPolyDataWriter(Input=vestecSamplingAlgorithm)
      coprocessor.RegisterWriter(xMLPPolyDataWriter1, filename='SpaceWeather_%t.pvtp', freq=1, paddingamount=0)
      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(appendDatasets1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

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
