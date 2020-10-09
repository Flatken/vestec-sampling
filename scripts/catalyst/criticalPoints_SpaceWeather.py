
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
# The following loads VESTEC's plugins.
# --------------------------------------------------------------
import os

## the path in which the Vestec Plugin is installed is different between linux and windows
## not clear why this is happening but we have to do the following workaround to correctly link
## the library files

if os.name == "posix":
	path_parent = os.path.dirname(os.getcwd())
	os.chdir(path_parent)
	vestec_plugin_path = os.path.join(path_parent,'lib/paraview-5.8/plugins/VestecPlugins/VestecPlugins.so')
else:	# windows
	vestec_plugin_path = os.path.join(os.getcwd(),'paraview-5.8/plugins/VestecPlugins/VestecPlugins.dll')
LoadPlugin(vestec_plugin_path, ns=globals())
print("Loading Plugin: "+vestec_plugin_path)


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

      # create a new 'VestecCriticalPointExtractionAlgorithm'
      vestecCriticalPointExtractionAlgorithm1 = VestecCriticalPointExtractionAlgorithm(Input=grid_)
      vestecCriticalPointExtractionAlgorithm1.Array = ['POINTS', 'B'] # for space-weather use-case

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(vestecCriticalPointExtractionAlgorithm1)
      # ----------------------------------------------------------------
	  
      # Now any catalyst writers
      writer = servermanager.writers.XMLPUnstructuredGridWriter(Input=vestecCriticalPointExtractionAlgorithm1)
      coprocessor.RegisterWriter(writer, filename='spaceWeather/VestecCriticalPointExtractionAlgorithm_%t.pvtu', freq=1, paddingamount=0)
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
