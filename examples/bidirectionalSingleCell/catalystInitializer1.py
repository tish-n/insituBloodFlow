from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.4.1 64 bits

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=2
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.4.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a producer from a simulation input
      cellspvd = coprocessor.CreateProducer(datadescription, 'cells')
      fluidpvd = coprocessor.CreateProducer(datadescription, 'fluid')
      dcpvd = coprocessor.CreateProducer(datadescription, 'dataCollection')
      # fluidpvd.CellArrays = ['vtkTestType']
      # fluidpvd.PointArrays = ['velocity', 'vorticity', 'velocityNorm', 'coords']

      # ----------------------------------------------------------------
      # finally, restore active source
      #   SetActiveSource(meshpvd)
      SetActiveSource(fluidpvd)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  # freqs = {'cells': [1],'fluid': [1], 'CLOTS':[1]}
  freqs = {'cells': [1],'fluid': [1],'dataCollection':[1]}
  # freqs = {'fluid': [1]}
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
coprocessor.EnableLiveVisualization(True, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    # datadescription.ForceOutput = True
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        print("GetForceOutput() is TRUE")
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
            print("SOMETHING has been added to datadescription")
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

    coprocessor.Pipeline.cellspvd.UpdatePipeline(datadescription.GetTime())
    coprocessor.Pipeline.fluidpvd.UpdatePipeline(datadescription.GetTime())
    coprocessor.Pipeline.dcpvd.UpdatePipeline(datadescription.GetTime())
    # coprocessor.Pipeline.clotpvd.UpdatePipeline(datadescription.GetTime())

    # coprocessor.Pipeline.cellspvd.UpdatePipeline()
    # coprocessor.Pipeline.fluidpvd.UpdatePipeline()

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
