import os
import vtk
import slicer
import qt
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from time import sleep
import numpy as np

#
# US_Cavity_Reconstruction
# collect multiple points for experiment point and create sphere
# adjust opacity -

class US_Cavity_Reconstruction(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "US Cavity Reconstruction"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["Ultrasound"]  # TODO: set categories (folders where the module shows up in the module selector)
    self.parent.dependencies = []  # TODO: add here list of module names that this module requires
    self.parent.contributors = ["Olivia Radcliffe (Perk/Medi Lab)"]  # TODO: replace with "Firstname Lastname (Organization)"
    # TODO: update with short description of the module and a link to online module documentation
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
See more information in <a href="https://github.com/organization/projectname#US_Cavity_Reconstruction">module documentation</a>.
"""
    # TODO: replace with organization, grant and thanks
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", registerSampleData)


#
# Register sample data sets in Sample Data module
#
# print (this is test data)
def registerSampleData():
  """
  Add data sets to Sample Data module.
  """
  # It is always recommended to provide sample data for users to make it easy to try the module,
  # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

  import SampleData
  iconsPath = os.path.join(os.path.dirname(__file__), 'Resources/Icons')

  # To ensure that the source code repository remains small (can be downloaded and installed quickly)
  # it is recommended to store data sets that are larger than a few MB in a Github release.

  # US_Cavity_Reconstruction1
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='US_Cavity_Reconstruction',
    sampleName='US_Cavity_Reconstruction1',
    # Thumbnail should have size of approximately 260x280 pixels and stored in Resources/Icons folder.
    # It can be created by Screen Capture module, "Capture all views" option enabled, "Number of images" set to "Single".
    thumbnailFileName=os.path.join(iconsPath, 'US_Cavity_Reconstruction1.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95",
    fileNames='US_Cavity_Reconstruction1.nrrd',
    # Checksum to ensure file integrity. Can be computed by this command:
    #  import hashlib; print(hashlib.sha256(open(filename, "rb").read()).hexdigest())
    checksums = 'SHA256:998cb522173839c78657f4bc0ea907cea09fd04e44601f17c82ea27927937b95',
    # This node name will be used when the data set is loaded
    nodeNames='US_Cavity_Reconstruction1'
  )

  # US_Cavity_Reconstruction2
  SampleData.SampleDataLogic.registerCustomSampleDataSource(
    # Category and sample name displayed in Sample Data module
    category='US_Cavity_Reconstruction',
    sampleName='US_Cavity_Reconstruction2',
    thumbnailFileName=os.path.join(iconsPath, 'US_Cavity_Reconstruction2.png'),
    # Download URL and target file name
    uris="https://github.com/Slicer/SlicerTestingData/releases/download/SHA256/1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97",
    fileNames='US_Cavity_Reconstruction2.nrrd',
    checksums = 'SHA256:1a64f3f422eb3d1c9b093d1a18da354b13bcf307907c66317e2463ee530b7a97',
    # This node name will be used when the data set is loaded
    nodeNames='US_Cavity_Reconstruction2'
  )


#
# US_Cavity_ReconstructionWidget
#

class US_Cavity_ReconstructionWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer).
    # Additional widgets can be instantiated manually and added to self.layout.
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/US_Cavity_Reconstruction.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = US_Cavity_ReconstructionLogic()

    # Connections


    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.tipSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.timerSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    
    
    # auto update checkbox
    self.observedTransformNode = None
    self.transformObserverTag = None
    self.ui.autoUpdateCheckBox.connect("toggled(bool)", self.onEnableAutoUpdate)
    self.ui.autoUpdateCheckBox.hide()

    # Buttons
    self.ui.loadModelsButton.connect('clicked(bool)', self.onloadModelsButton)
    self.ui.placeTestPointButton.connect('clicked(bool)', self.onPlaceTestPointButton)
    self.ui.placeNeedleButton.connect('clicked(bool)', self.onPlaceNeedleButton)
    self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.ui.stopButton.connect('clicked(bool)', self.onStopButton)
    self.ui.generateModel.connect('clicked(bool)', self.onGenerateButton)

    # visualization
    self.ui.generateTumourButton.connect('clicked(bool)', self.onGenerateTumorButton)
    self.ui.enableNavigationCheckbox.connect("toggled(bool)", self.onEnableNavigationCheckbox)


    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()

  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()
  
  def updateGUIFromParameterNode(self, caller=None, event=None):
    if self._parameterNode is None:
      return
    self.ui.applyButton.enabled = self._parameterNode.GetNodeReference("OutputPoints")
    self.ui.generateModel.enabled = self._parameterNode.GetNodeReference("OutputPoints")

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    # Update node selectors and sliders
    self.ui.tipSelector.setCurrentNode(self._parameterNode.GetNodeReference("ProbeTip"))
    self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputPoints"))
    #self.ui.timerSelector.setCurrentNode(self._parameterNode.GetNodeReference("TimerSeconds"))
    
    # Update buttons states and tooltips
    if self._parameterNode.GetNodeReference("ProbeTip") and self._parameterNode.GetNodeReference("OutputPoints"):
      self.ui.applyButton.toolTip = "Compute output model"
      self.ui.applyButton.enabled = True
    else:
      self.ui.applyButton.toolTip = "Select input and output point lists"
      self.ui.applyButton.enabled = False

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID("ProbeTip", self.ui.tipSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID("OutputPoints", self.ui.outputSelector.currentNodeID)
    #self._parameterNode.SetNodeReferenceID("TimerSeconds", self.ui.timerSelector.currentNodeID)
    

    self._parameterNode.EndModify(wasModified)
    
  def onEnableAutoUpdate(self, autoUpdate):
    if self.transformObserverTag:
      self.observedTransformNode.RemoveObserver(self.transformObserverTag)
      self.observedTransformNode = None
      self.transformObserverTag = None

    if autoUpdate: 
      self.observedTransformNode = slicer.util.getNode("ProbeToTracker")
      self.transformObserverTag = self.observedTransformNode.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, self.onTransformNodeModified)
  
  def onTransformNodeModified(self, caller=None, event=None):
    self.collectedPoints = self.logic.placePoint(self.ui.tipSelector.currentNode(), self.ui.outputSelector.currentNode())

  def onloadModelsButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
      probe = slicer.util.loadModel(os.path.dirname(__file__) +"\Resources\\ExampleModels/new_probe.stl")
      retractor = slicer.util.loadModel(os.path.dirname(__file__) +"\Resources\\ExampleModels/new_retractor.stl")
  
  def onPlaceTestPointButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
      import numpy as np

      probeTip = self.ui.tipSelector.currentNode()

      pos = np.zeros(3)
      probeTip.GetNthFiducialPosition(0,pos)

      tip = slicer.mrmlScene.GetFirstNodeByName("tip")
      if tip is None:
        tip = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "tip")
      
      tip.AddControlPoint(pos)

      tip.SetDisplayVisibility(0)

      testPoints = slicer.mrmlScene.GetFirstNodeByName("testPoints")
      if testPoints is None:
        testPoints = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", "testPoints")
      testPoints.GetDisplayNode().SetSelectedColor(0.3, 0.6, 0.1)
      testPoints.GetDisplayNode().UseGlyphScaleOff()
      testPoints.GetDisplayNode().SetGlyphSize(5)

      RetractorToTrackerNode = slicer.util.getNode("RetractorToTracker")
      testPoints.SetAndObserveTransformNodeID(RetractorToTrackerNode.GetID())

      self.logic.placePoint(tip, testPoints)

      print(testPoints.GetNumberOfMarkups() )

            
      for i in range(testPoints.GetNumberOfMarkups()):

        import numpy as np
        pos =  np.zeros(3)
        testPoints.GetNthFiducialPosition(0,pos)

        pointModel = self.logic.generateTestPointModel(pos, "testPointModel")

        transformID = slicer.util.getNode("NeedleTipToNeedle").GetID()
        self.logic.addBreachWarning(pointModel, transformID)
  
  def onPlaceNeedleButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
      needleModel = slicer.mrmlScene.GetFirstNodeByName("NeedleModel")
      transform = needleModel.GetParentTransformNode()

      insertedNeedleModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode")
      insertedNeedleModel.CopyContent(needleModel)
      insertedNeedleModel.CreateDefaultDisplayNodes()
      insertedNeedleModel.SetDisplayVisibility(1)

      insertedNeedleModel.SetAndObserveTransformNodeID(transform.GetID())

      insertedNeedleModel.HardenTransform()
      
  
  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      time = self.ui.timerSelector.value

      """self.timer = qt.QTimer()
      self.timer.timeout.connect(self.timeout())
      #self.timer.timeout.connect(qt.processOneThing) 

      self.timer.start(time*1000)

"""
      #self.ui.countdownValueLabel.text = str(self.timer.remainingTime())
      self.ui.countdownValueLabel.text = str(time)
      self.ui.applyButton.enabled = False
      self.ui.stopButton.enabled = True
      
      # TODO: fix so it shows in UI!
      while time:
        sleep(1)
        time -= 1
        self.ui.countdownValueLabel.text = str(time)
      
      RetractorToTrackerNode = slicer.util.getNode("RetractorToTracker")
      self.ui.outputSelector.currentNode().SetAndObserveTransformNodeID(RetractorToTrackerNode.GetID())

      self.ui.autoUpdateCheckBox.checked = True
      self.collectedPoints = self.logic.placePoint(self.ui.tipSelector.currentNode(), self.ui.outputSelector.currentNode())

  """def timeout(self):
    if self.timer.remainingTime == 0:
      self.timer.stop
    self.ui.countdownValueLabel.text = str(self.timer.remainingTime())
  """

  def onStopButton(self):
    """
    Stop running processing when user clicks "Stop" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      #stop collecting
      self.ui.autoUpdateCheckBox.checked = False
      self.ui.applyButton.enabled = True

  def onGenerateButton(self):
    
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      self.logic.process(self.ui.tipSelector.currentNode(), 
      self.ui.outputSelector.currentNode())

  def onGenerateTumorButton(self):
    """
    Generates a tumor model from the pointListRed_World points
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
      # TODO: add selector so not hardcoded?
      self.logic.generateTumor(slicer.util.getNode("pointListRed_World"))

  
  def onEnableNavigationCheckbox(self, state):
    """
    Enable/Disable navigation
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
      self.logic.enableNavigation(state)

#
# US_Cavity_ReconstructionLogic
#

class US_Cavity_ReconstructionLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)

    self.num = 0

  def setDefaultParameters(self, parameterNode):
    """
    Initialize parameter node with default settings.
    """
    if not parameterNode.GetParameter("Threshold"):
      parameterNode.SetParameter("Threshold", "100.0")
    if not parameterNode.GetParameter("Invert"):
      parameterNode.SetParameter("Invert", "false")
  
  def placePoint(self, probeTip, outputPoints):
    import numpy as np
    print("Debugging: point placed")
    for i in range(probeTip.GetNumberOfMarkups()):
      # get point position
      pos = np.zeros(3)
      probeTip.GetNthFiducialPosition(i,pos)

      #TODO: get transforms that the model are under rather than getting them by name?
      # get and apply probe transforms to point position in retractor coordinates
      probeModelToProbeNode = slicer.util.getNode("ProbeModelToProbe")
      probeModelToProbeMatrix = probeModelToProbeNode.GetMatrixTransformToParent()
      pos = probeModelToProbeMatrix.MultiplyPoint(np.append(pos,1))

      ProbeToRetractorNode = slicer.util.getNode("ProbeToRetractor")
      ProbeToRetractorMatrix = ProbeToRetractorNode.GetMatrixTransformToParent()
      pos = ProbeToRetractorMatrix.MultiplyPoint(pos)

      # add point to output point list
      n = outputPoints.AddControlPoint(pos[:3])
      #outputPoints.SetNthControlPointLabel(n, str(self.num))
      outputPoints.SetNthControlPointLabel(n, "")
      # set the visibility flag
      outputPoints.SetNthControlPointVisibility(n, 1)

      return outputPoints


  def generateModel(self, pointCloud, name):

    import numpy as np

    pointsForHull = vtk.vtkPoints()
    Zs = []

    # add z value to zs and positions to hull points
    for i in range(pointCloud.GetNumberOfMarkups()):
        pos = np.zeros(3)
        pointCloud.GetNthFiducialPosition(i,pos)

        pointsForHull.InsertNextPoint(pos)
        Zs.append(pos[2])

    # average the 5 highest points to determine heihght of the phantom
    Zs.sort(reverse=True)
    sum = 0
    for val in Zs[:50]:
      sum += val

    height = sum/50

    print(height)

    phantomHeight = 50

    # create convex hull of the cavity
    hullPolydata = vtk.vtkPolyData()
    hullPolydata.SetPoints(pointsForHull)

    delaunay = vtk.vtkDelaunay3D()
    delaunay.SetInputData(hullPolydata)
    delaunay.Update()

    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputConnection(delaunay.GetOutputPort())
    surfaceFilter.Update()

    subdivisionFilter = vtk.vtkButterflySubdivisionFilter()
    subdivisionFilter.SetInputConnection(surfaceFilter.GetOutputPort())
    subdivisionFilter.SetNumberOfSubdivisions(3)
    subdivisionFilter.Update()

    convexHull = vtk.vtkDelaunay3D()
    convexHull.SetInputConnection(subdivisionFilter.GetOutputPort())
    convexHull.Update()

    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputData(convexHull.GetOutput())
    surfaceFilter.Update()

    # create convex hull model
    outputModel = slicer.mrmlScene.GetFirstNodeByName(name)
    if outputModel is None:
      outputModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", name)
    outputModel.SetAndObservePolyData(surfaceFilter.GetOutput())
    outputModel.CreateDefaultDisplayNodes()
    outputModel.GetDisplayNode().SetColor(0,0,1)

    #TODO: get transform retractor is under rathe rthan by name
    # transform cavity model to retractor space so that it moves with the retractor
    retractorToTrackerNode = slicer.util.getNode("RetractorToTracker")
    outputModel.SetAndObserveTransformNodeID(retractorToTrackerNode.GetID())

   
    return outputModel
    

  def setupForSubtract(self, tumorModel):
    """
    Create all elements for subtraction using the segmentation editor modeule
    :param breastModel: output model - model an object is subtracted from
    :param tumorModel: model to subtract
    :param segmentationNode: an optional pre-existing segmentation Node
    """
    # TODO: adjust way of getting nodes to be more universal

     # load in the breast pphantom model and transform it to the retractor coordinate space and proper height 
    breastModel = slicer.util.loadModel(os.path.dirname(__file__) +"\Resources\\ExampleModels/breastPhantom.stl")

    #breastModel = slicer.util.getNode("breastPhantom_1")
    cavityTopToRetractorNode = slicer.mrmlScene.GetFirstNodeByName("CavityTopToRetractor")
    if cavityTopToRetractorNode is None:
      breastPhantomToRetractorFilename = os.path.dirname(__file__) + "\Resources\\CavityTopToRetractor_new_2.h5"
      cavityTopToRetractorNode = slicer.util.loadTransform(breastPhantomToRetractorFilename)

    """
    breastPhantomToCavityTopNode = slicer.mrmlScene.GetFirstNodeByName("BreastPhantomToCavityTop")
    if breastPhantomToCavityTopNode is None:
      breastPhantomToCavityTopNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode", "BreastPhantomToCavityTop")

    breastPhantomToCavityTopTransform = vtk.vtkTransform()
    breastPhantomToCavityTopTransform.Translate(0, 0, (height+(phantomHeight/2)))
    breastPhantomToCavityTopNode.SetMatrixTransformToParent(breastPhantomToCavityTopTransform.GetMatrix())"""

    retractorToTrackerNode = slicer.util.getNode("RetractorToTracker")

    breastModel.SetAndObserveTransformNodeID(cavityTopToRetractorNode.GetID())
    #breastPhantomToCavityTopNode.SetAndObserveTransformNodeID(cavityTopToRetractorNode.GetID())
    cavityTopToRetractorNode.SetAndObserveTransformNodeID(retractorToTrackerNode.GetID())


    breastModel.SetDisplayVisibility(0)
    
    segmentationNode = slicer.mrmlScene.GetFirstNodeByName("CavitySegmentation")
    if segmentationNode is None:
      segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "CavitySegmentation")
    segmentationNode.CreateDefaultDisplayNodes()
    
    #TODO get retractor transform rather than using name
    retractorModelToRetractorNode = slicer.util.getNode("RetractorModelToRetractor")
    segmentationNode.SetAndObserveTransformNodeID(retractorModelToRetractorNode.GetID())

    # import models to segmentations
    slicer.modules.segmentations.logic().ImportModelToSegmentationNode(breastModel, segmentationNode)

    #tumor = slicer.util.getNode("solidCavity")
    tumorModel.SetDisplayVisibility(0)

    slicer.modules.segmentations.logic().ImportModelToSegmentationNode(tumorModel, segmentationNode)

    # create segment editor nodes
    segmentEditorNode = slicer.vtkMRMLSegmentEditorNode()
    slicer.mrmlScene.AddNode(segmentEditorNode)

    segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
    segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
    segmentEditorWidget.setSegmentationNode(segmentationNode)
    segmentEditorWidget.setMRMLScene(slicer.mrmlScene)

    breastID = segmentationNode.GetSegmentation().GetSegmentIdBySegmentName(breastModel.GetName())
    tumorID = segmentationNode.GetSegmentation().GetSegmentIdBySegmentName(tumorModel.GetName())

    slicer.mrmlScene.RemoveNode(breastModel)

    return segmentEditorWidget, segmentEditorNode, segmentationNode, breastID, tumorID


  def process(self, probeTip, outputPoints, showResult=True):
    """
    Run the processing algorithm.s
    Can be used without GUI widget.
    :param probeTip: points on probe tip used for scan
    :param outputPoints: scanned point cloud
    :param showResult: show output volume in slice viewers
    """
    #self.pointCloud = self.generatePointCloud(outputPoints)
    # hide users point cloud
    outputPoints.SetDisplayVisibility(0)
    tumorModel = self.generateModel(outputPoints, "solidCavity")
    segmentEditorWidget, segmentEditorNode, segmentationNode, breastID, tumorID = self.setupForSubtract(tumorModel)
    self.subtract_segment(segmentEditorWidget, segmentEditorNode, segmentationNode, breastID, tumorID)
    self.set_triple3D()

  def subtract_segment(self, segmentEditorWidget, segmentEditorNode, segmentationNode, segmentID, to_subtract_segmentID):
    '''Perform logical subtraction of to_subtract_segmentID from segmentID'''
    import SegmentEditorEffects
    
    segmentEditorWidget.setActiveEffectByName('Logical operators')

    effect = segmentEditorWidget.activeEffect()
    effect.setParameter("Operation", SegmentEditorEffects.LOGICAL_SUBTRACT)
    effect.setParameter("BypassMasking",1)
    effect.setParameter("ModifierSegmentID",to_subtract_segmentID)
    effect.self().onApply()

    segmentationNode.GetDisplayNode().SetSegmentVisibility(to_subtract_segmentID, 0)
    segmentationNode.RemoveSegment(to_subtract_segmentID)


    slicer.mrmlScene.RemoveNode(segmentEditorNode)

    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    exportFolderItemId = shNode.CreateFolderItem(shNode.GetSceneItemID(), "CavityVisualization")
  
    slicer.modules.segmentations.logic().ExportAllSegmentsToModels(segmentationNode, exportFolderItemId)
    
    cavity = slicer.mrmlScene.GetFirstNodeByName((segmentationNode.GetSegmentation().GetSegment(segmentID)).GetName())
    cavity.GetDisplayNode().SetOpacity(0.5)

    self.removeSegmentation(segmentationNode, segmentID)



  def removeSegmentation(self, segmentationNode, segmentID):
    segmentationDisplay = segmentationNode.GetDisplayNode()
    segmentationDisplay.SetSegmentVisibility(segmentID, 0)
    
    slicer.mrmlScene.RemoveNode(segmentationDisplay)
    slicer.mrmlScene.RemoveNode(segmentationNode)

  def set_triple3D(self):
    # Set to triple 3D view
    triple3D = slicer.vtkMRMLLayoutNode.SlicerLayoutTriple3DEndoscopyView
    slicer.app.layoutManager().setLayout(triple3D)
    # Set the angle of the 3D views
    
    # The second view should be the side view
    # The third view should be the front view
    cameraNode = slicer.modules.cameras.logic().GetViewActiveCameraNode(slicer.app.layoutManager().threeDWidget(0).mrmlViewNode())
    # Set the camera position 1
    cameraNode.SetPosition([-22.797165793169604, -327.31100424393105, -32.89968700041658])
    cameraNode.SetFocalPoint([-6.428269131321022, 25.387306384949536, -22.541016228068706])
    cameraNode.SetViewUp([-1, 0, 0])
    # Set the camera position 2
    cameraNode = slicer.modules.cameras.logic().GetViewActiveCameraNode(slicer.app.layoutManager().threeDWidget(1).mrmlViewNode())
    cameraNode.SetPosition([-29.739841516419293, 6.800195873066253, -281.50844451170815])
    cameraNode.SetFocalPoint([0.5996254445084617, 22.754982645586885, -1.3608656337910645])
    cameraNode.SetViewUp([0, -1, 0])
    # Set the camera position 3
    cameraNode = slicer.modules.cameras.logic().GetViewActiveCameraNode(slicer.app.layoutManager().threeDWidget(2).mrmlViewNode())
    cameraNode.SetPosition([224.2681082329677, 14.389053723089877, -75.35438455128524])
    cameraNode.SetFocalPoint([-8.528836308582385, 17.923943384488755, -61.19913412126311])
    cameraNode.SetViewUp([-0.010467168505081623, -0.9969893997410547, 0.0768282186924673])

  def generateTestPointModel(self, center, name):
  
    outputModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", name)
    
    s = vtk.vtkSphereSource()
    s.SetRadius(2.5)
    s.SetCenter(center)
    s.Update()

    outputModel.SetAndObservePolyData(s.GetOutput())
    outputModel.CreateDefaultDisplayNodes()
    outputModel.GetDisplayNode().SetColor(0,0,1)

    RetractorToTrackerNode = slicer.util.getNode("RetractorToTracker")
    outputModel.SetAndObserveTransformNodeID(RetractorToTrackerNode.GetID())

    return outputModel

  def addBreachWarning(self, pointModel, transformID):
    breachWarningNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLBreachWarningNode")
    slicer.modules.breachwarning.logic().SetWatchedModelNode(pointModel,breachWarningNode)
    breachWarningNode.SetAndObserveToolTransformNodeId(transformID)
    slicer.modules.breachwarning.logic().SetLineToClosestPointVisibility(1,breachWarningNode)

  def generateTumor(self, pointCloud):
    pointsForHull = vtk.vtkPoints()
    Zs = []

    # add z value to zs and positions to hull points
    for i in range(pointCloud.GetNumberOfMarkups()):
        pos = np.zeros(3)
        pointCloud.GetNthFiducialPosition(i,pos)

        pointsForHull.InsertNextPoint(pos)
        Zs.append(pos[2])

    pointCloud.SetDisplayVisibility(0)

    # create convex hull of the tumor
    hullPolydata = vtk.vtkPolyData()
    hullPolydata.SetPoints(pointsForHull)

    delaunay = vtk.vtkDelaunay3D()
    delaunay.SetInputData(hullPolydata)
    delaunay.Update()

    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputConnection(delaunay.GetOutputPort())
    surfaceFilter.Update()

    subdivisionFilter = vtk.vtkButterflySubdivisionFilter()
    subdivisionFilter.SetInputConnection(surfaceFilter.GetOutputPort())
    subdivisionFilter.SetNumberOfSubdivisions(3)
    subdivisionFilter.Update()

    convexHull = vtk.vtkDelaunay3D()
    convexHull.SetInputConnection(subdivisionFilter.GetOutputPort())
    convexHull.Update()

    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputData(convexHull.GetOutput())
    surfaceFilter.Update()

    # create convex hull model
    outputModel = slicer.mrmlScene.GetFirstNodeByName("Tumour")
    if outputModel is None:
      outputModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", "Tumour")
    outputModel.SetAndObservePolyData(surfaceFilter.GetOutput())
    outputModel.CreateDefaultDisplayNodes()
    outputModel.GetDisplayNode().SetColor(1,0,0)
  
  def enableNavigation(self, state):
    if state:
      tumorModel = slicer.util.getNode("Tumour") 

      transformID = slicer.util.getNode("ProbeModelToProbe").GetID()
      self.addBreachWarning(tumorModel, transformID)
    else:
      breachWarningNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLBreachWarningNode")
      slicer.modules.breachwarning.logic().SetWatchedModelNode(None,breachWarningNode)

      distance = slicer.util.getNode("d")
      slicer.mrmlScene.RemoveNode(distance)


    


#
# US_Cavity_ReconstructionTest
#

class US_Cavity_ReconstructionTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_US_Cavity_Reconstruction1()

  def test_US_Cavity_Reconstruction1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    # Get/create input data

    import SampleData
    registerSampleData()
    inputModel = SampleData.downloadSample('US_Cavity_Reconstruction1')
    self.delayDisplay('Loaded test data set')

    inputScalarRange = inputModel.GetImageData().GetScalarRange()
    self.assertEqual(inputScalarRange[0], 0)
    self.assertEqual(inputScalarRange[1], 695)

    probeTip = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    outputPoints = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")

    # Test the module logic

    logic = US_Cavity_ReconstructionLogic()

    # Test algorithm with non-inverted threshold
    logic.process(outputPoints, True)
    outputScalarRange = outputPoints.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(outputPoints, False)
    outputScalarRange = outputPoints.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
