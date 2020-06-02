
#--------------------------------------------------------------

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

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.6.3
#--------------------------------------------------------------

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

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.6.3

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.6.3
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      grid_ = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Tetrahedralize'
      tetrahedralize1 = Tetrahedralize(Input=grid_)

      # create a new 'VestecCriticalPointExtractionAlgorithm'
      vestecCriticalPointExtractionAlgorithm1 = VestecCriticalPointExtractionAlgorithm(Input=tetrahedralize1)
      vestecCriticalPointExtractionAlgorithm1.Array = ['POINTS', 'Vec']

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(vestecCriticalPointExtractionAlgorithm1)
      # ----------------------------------------------------------------
	  
      # Now any catalyst writers
      xMLPPolyDataWriter1 = servermanager.writers.XMLPPolyDataWriter(Input=vestecCriticalPointExtractionAlgorithm1)
      coprocessor.RegisterWriter(xMLPPolyDataWriter1, filename='results/VestecCriticalPointExtractionAlgorithm_%t.pvtp', freq=1, paddingamount=0)
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
