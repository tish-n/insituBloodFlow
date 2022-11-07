# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
from paraview import coprocessing
# #### disable automatic camera reset on 'Show'
# paraview.simple._DisableFirstRenderCameraReset()

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
      #   oscillators = coprocessor.CreateProducer(datadescription, 'oscillators')
      #   meshpvd = coprocessor.CreateProducer(datadescription, 'mesh')
      cellspvd = coprocessor.CreateProducer(datadescription, 'cells')  
      fluidpvd = coprocessor.CreateProducer(datadescription, 'fluid')
      # ----------------------------------------------------------------
      # finally, restore active source
      # SetActiveSource(meshpvd)
      SetActiveSource(cellspvd)
      SetActiveSource(fluidpvd)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  # freqs = {'oscillators': [1], 'mesh': [1], 'fluid'}
  freqs = {'fluid': [1], 'cells': [1]}
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

    coprocessor.Pipeline.meshpvd.UpdatePipeline(datadescription.GetTime())
    coprocessor.Pipeline.oscillators.UpdatePipeline(datadescription.GetTime())

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

# remainder

# # Create a new 'Render View'
# renderView1 = CreateView('RenderView')
# renderView1.ViewSize = [2267, 1173]
# renderView1.AxesGrid = 'GridAxes3DActor'
# renderView1.CenterOfRotation = [15.0, 15.0, 40.0]
# renderView1.StereoType = 0
# renderView1.CameraPosition = [140.49911592578337, 24.493254706733072, 79.36051468750458]
# renderView1.CameraFocalPoint = [15.000000000000016, 14.999999999999966, 40.00000000000004]
# renderView1.CameraViewUp = [0.006218855007531415, 0.9674058754836297, -0.25315449417726577]
# renderView1.CameraParallelScale = 49.969990994595946
# renderView1.Background = [0.32, 0.34, 0.43]

# # ----------------------------------------------------------------
# # setup the data processing pipelines
# # ----------------------------------------------------------------

# # create a new 'PVD Reader'
# cellspvd = PVDReader(FileName='/home/tishn/myFork/singleCell4.0/vtk_output/cells')

# # create a new 'PVD Reader'
# fluidpvd = PVDReader(FileName='/home/tishn/myFork/singleCell4.0/vtk_output/fluid')

# # ----------------------------------------------------------------
# # setup color maps and opacity mapes used in the visualization
# # note: the Get..() functions create a new object, if needed
# # ----------------------------------------------------------------

# # get color transfer function/color map for 'vtkBlockColors'
# vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
# vtkBlockColorsLUT.InterpretValuesAsCategories = 1
# vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
# vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1', '2', '3']
# vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# # get opacity transfer function/opacity map for 'vtkBlockColors'
# vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# # ----------------------------------------------------------------
# # setup the visualization in view 'renderView1'
# # ----------------------------------------------------------------

# # show data from fluidpvd
# fluidpvdDisplay = Show(fluidpvd, renderView1)
# # trace defaults for the display properties.
# fluidpvdDisplay.Representation = 'Outline'
# fluidpvdDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
# fluidpvdDisplay.LookupTable = vtkBlockColorsLUT
# fluidpvdDisplay.OSPRayScaleArray = 'velocity'
# fluidpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
# fluidpvdDisplay.SelectOrientationVectors = 'None'
# fluidpvdDisplay.ScaleFactor = 8.6
# fluidpvdDisplay.SelectScaleArray = 'None'
# fluidpvdDisplay.GlyphType = 'Arrow'
# fluidpvdDisplay.GlyphTableIndexArray = 'None'
# fluidpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
# fluidpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
# fluidpvdDisplay.GaussianRadius = 4.3
# fluidpvdDisplay.SetScaleArray = ['POINTS', 'velocityNorm']
# fluidpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
# fluidpvdDisplay.OpacityArray = ['POINTS', 'velocityNorm']
# fluidpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# # show color legend
# fluidpvdDisplay.SetScalarBarVisibility(renderView1, True)

# # show data from cellspvd
# cellspvdDisplay = Show(cellspvd, renderView1)
# # trace defaults for the display properties.
# cellspvdDisplay.Representation = 'Surface'
# cellspvdDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
# cellspvdDisplay.LookupTable = vtkBlockColorsLUT
# cellspvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
# cellspvdDisplay.SelectOrientationVectors = 'None'
# cellspvdDisplay.ScaleFactor = 1.7999815078354269
# cellspvdDisplay.SelectScaleArray = 'None'
# cellspvdDisplay.GlyphType = 'Arrow'
# cellspvdDisplay.GlyphTableIndexArray = 'None'
# cellspvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
# cellspvdDisplay.PolarAxes = 'PolarAxesRepresentation'
# cellspvdDisplay.GaussianRadius = 0.8999907539177134
# cellspvdDisplay.SetScaleArray = [None, '']
# cellspvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
# cellspvdDisplay.OpacityArray = [None, '']
# cellspvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# # show color legend
# cellspvdDisplay.SetScalarBarVisibility(renderView1, True)

# # setup the color legend parameters for each legend in this view

# # get color legend/bar for vtkBlockColorsLUT in view renderView1
# vtkBlockColorsLUTColorBar = GetScalarBar(vtkBlockColorsLUT, renderView1)
# vtkBlockColorsLUTColorBar.Title = 'vtkBlockColors'
# vtkBlockColorsLUTColorBar.ComponentTitle = ''

# # ----------------------------------------------------------------
# # finally, restore active source
# SetActiveSource(fluidpvd)
# # ----------------------------------------------------------------