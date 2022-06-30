import os
import vtk
import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
from time import sleep

#
# US_Cavity_Reconstruction
#

class US_Cavity_Reconstruction(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "US_Cavity_Reconstruction"  # TODO: make this more human readable by adding spaces
    self.parent.categories = ["Examples"]  # TODO: set categories (folders where the module shows up in the module selector)
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

    #self.MarkupsToModelLogic=slicer.modules.markupstomodel.logic()

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.modelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
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
    self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.ui.stopButton.connect('clicked(bool)', self.onStopButton)
    self.ui.generateModel.connect('clicked(bool)', self.onGenerateButton)

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

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    if not self._parameterNode.GetNodeReference("InputModel"):
      firstModelNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLModelNode")
      if firstModelNode:
        self._parameterNode.SetNodeReferenceID("InputModel", firstModelNode.GetID())
    

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
    self.ui.applyButton.enabled = self._parameterNode.GetNodeReference("InputModel")
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
    self.ui.modelSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputModel"))
    self.ui.tipSelector.setCurrentNode(self._parameterNode.GetNodeReference("ProbeTip"))
    self.ui.outputSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputPoints"))
    #self.ui.timerSelector.setCurrentNode(self._parameterNode.GetNodeReference("TimerSeconds"))
    
    # Update buttons states and tooltips
    if self._parameterNode.GetNodeReference("InputModel") and self._parameterNode.GetNodeReference("OutputPoints"):
      self.ui.applyButton.toolTip = "Compute output volume"
      self.ui.applyButton.enabled = True
    else:
      self.ui.applyButton.toolTip = "Select input and output volume nodes"
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

    self._parameterNode.SetNodeReferenceID("InputModel", self.ui.modelSelector.currentNodeID)
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
      probe = slicer.util.loadModel(os.path.dirname(__file__) +"\Resources\\ExampleModels/built_probe.stl")
      retractor = slicer.util.loadModel(os.path.dirname(__file__) +"\Resources\\ExampleModels/retractorStand_v3.stl")
  
  
  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      t = self.ui.timerSelector.value

      self.ui.countdownValueLabel.text = str(t)
      self.ui.applyButton.enabled = False
      self.ui.stopButton.enabled = True
      
      # TODO: fix so it shows in UI!
      while t:
        sleep(1)
        t -= 1
        self.ui.countdownValueLabel.text = str(t)
      
      t = slicer.util.getNode("RetractorToTracker")
      self.ui.outputSelector.currentNode().SetAndObserveTransformNodeID(t.GetID())

      self.ui.autoUpdateCheckBox.checked = True
      self.collectedPoints = self.logic.placePoint(self.ui.tipSelector.currentNode(), self.ui.outputSelector.currentNode())
  
  def onStopButton(self):
    """
    Stop running processing when user clicks "Stop" button.
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      #stop collecting
      self.ui.autoUpdateCheckBox.checked = False
      self.ui.applyButton.enabled = True
      #self.removeObservers()
      print("should stop")

  def onGenerateButton(self):
    """
    .
    """
    with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):

      self.logic.process(self.ui.modelSelector.currentNode(),
      self.ui.tipSelector.currentNode(), 
      self.ui.outputSelector.currentNode())

      

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

    # TODO: transform points to retractor so that they move with the retractor
  
    self.num +=1

    for i in range(probeTip.GetNumberOfMarkups()):
      pos = np.zeros(3)
      probeTip.GetNthFiducialPosition(i,pos)

      transform = slicer.util.getNode("probeModelToProbe")
      matrix = transform.GetMatrixTransformToParent()
      pos = matrix.MultiplyPoint(np.append(pos,1))

      transform2 = slicer.util.getNode("ProbeToRetractor")
      matrix2 = transform2.GetMatrixTransformToParent()
      pos = matrix2.MultiplyPoint(pos)

      n = outputPoints.AddControlPoint(pos[:3])
      outputPoints.SetNthControlPointLabel(n, str(self.num))
      # set the visibility flag
      outputPoints.SetNthControlPointVisibility(n, 1)

      return outputPoints

    """pos = np.zeros(3)
    probeTip.GetNthFiducialPosition(0,pos)

    transform = slicer.util.getNode("probeModelToProbe")
    matrix = transform.GetMatrixTransformToParent()
    pos = matrix.MultiplyPoint(np.append(pos,1))

    transform2 = slicer.util.getNode("ProbeToTracker")
    matrix2 = transform2.GetMatrixTransformToParent()
    pos = matrix2.MultiplyPoint(pos)

    n = outputPoints.AddControlPoint(pos[:3])
    outputPoints.SetNthControlPointLabel(n, str(self.num))
    # set the visibility flag
    outputPoints.SetNthControlPointVisibility(n, 1)"""

  """def generatePointCloud(self,outputPoints):

      # hide collected points
      outputPoints.SetDisplayVisibility(0)

      import numpy as np

      # go through collected point to find the outermost points
      for i in range(outputPoints.GetNumberOfMarkups()):
        pos = np.zeros(3)
        outputPoints.GetNthFiducialPosition(i,pos)

        # initialize values to first point
        if i == 0:
          highZ = pos[2]
          lowZ = pos[2]
          highZPos = pos
          lowZPos = pos

          highY = pos[1]
          lowY = pos[1]
          highYPos = pos
          lowYPos = pos

          highX = pos[0]
          lowX = pos[0]
          highXPos = pos
          lowXPos = pos
        
        else:
          if pos[2] > highZ:
            highZ = pos[2]
            highZPos = pos
          if pos[2] < lowZ:
            lowZ = pos[2]
            lowZPos = pos
          
          if pos[1] > highY:
            highY = pos[1]
            highYPos = pos
          if pos[1] < lowY:
            lowY = pos[1]
            lowYPos = pos

          if pos[0] > highX:
            highX = pos[0]
            highXPos = pos
          if pos[0] < lowX:
            lowX = pos[0]
            lowXPos = pos

      # create new fuducial list and add/show outermost points
      pointListNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      # TODO: change the name of the fudical list
      pointListNode.SetName("Cavity")
      n1 = pointListNode.AddControlPoint(highZPos)
      n2 = pointListNode.AddControlPoint(lowZPos)
      pointListNode.SetNthControlPointVisibility(n1, 1)
      pointListNode.SetNthControlPointVisibility(n2, 1)

      n3 = pointListNode.AddControlPoint(highYPos)
      n4 = pointListNode.AddControlPoint(lowYPos)
      pointListNode.SetNthControlPointVisibility(n3, 1)
      pointListNode.SetNthControlPointVisibility(n4, 1)

      n5 = pointListNode.AddControlPoint(highXPos)
      n6 = pointListNode.AddControlPoint(lowXPos)
      pointListNode.SetNthControlPointVisibility(n5, 1)
      pointListNode.SetNthControlPointVisibility(n6, 1)

      return pointListNode"""

  def generateModel(self, pointList, inputModel):

    import numpy as np

    pointList.SetDisplayVisibility(0)

    pointsForHull = vtk.vtkPoints()
    
    for i in range(pointList.GetNumberOfMarkups()):
        pos = np.zeros(3)
        pointList.GetNthFiducialPosition(i,pos)
        pointsForHull.InsertNextPoint(pos)

        if i == 0:
          highestZ = pos[2]
        else:
          if pos[2] > highestZ:
            highestZ = pos[2]

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

    outputModel = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode")
    outputModel.SetAndObservePolyData(surfaceFilter.GetOutput())
    outputModel.SetName("solidCavity")
    outputModel.CreateDefaultDisplayNodes()
    outputModel.GetDisplayNode().SetColor(0,0,1)

    
    t = slicer.util.getNode("RetractorToTracker")
    outputModel.SetAndObserveTransformNodeID(t.GetID())

    breast = slicer.util.loadModel(os.path.dirname(__file__) +"\Resources\\ExampleModels/breastPhantom.stl")

    bTransform = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode")
    matrix = vtk.vtkMatrix4x4()
    matrix.SetElement(0, 0, -1.28)
    matrix.SetElement(1, 1, 1.0)
    matrix.SetElement(2, 2, -1.28)
    matrix.SetElement(0, 2, -0.1)
    matrix.SetElement(2, 0, 0.1)
    matrix.SetElement(0, 3, 11.35)
    matrix.SetElement(1, 3, 51.0)
    matrix.SetElement(2, 3, highestZ*2.4)
    bTransform.SetMatrixTransformToParent(matrix)

    breast.SetAndObserveTransformNodeID(bTransform.GetID())
    bTransform.SetAndObserveTransformNodeID(t.GetID())
  

    segmentEditorWidget, segmentEditorNode, segmentationNode, breastID, tumorID = self.setupForSubtract(breast)
    self.subtract_segment(segmentEditorWidget, segmentEditorNode, segmentationNode, breastID, tumorID)


  def setupForSubtract(self, inputModel, segmentationNode=None):
    # TODO: adjust way of getting nodes to be more universal

    inputModel.SetDisplayVisibility(0)
    
    segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
    segmentationNode.CreateDefaultDisplayNodes()
    
    transform = slicer.util.getNode("retractorModelToRetractor")
    segmentationNode.SetAndObserveTransformNodeID(transform.GetID())

    slicer.modules.segmentations.logic().ImportModelToSegmentationNode(inputModel, segmentationNode)

    tumor = slicer.util.getNode("solidCavity")
    tumor.SetDisplayVisibility(0)

    slicer.modules.segmentations.logic().ImportModelToSegmentationNode(tumor, segmentationNode)

    segmentEditorNode = slicer.vtkMRMLSegmentEditorNode()
    slicer.mrmlScene.AddNode(segmentEditorNode)

    segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
    segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
    segmentEditorWidget.setSegmentationNode(segmentationNode)
    segmentEditorWidget.setMRMLScene(slicer.mrmlScene)

    breastID = segmentationNode.GetSegmentation().GetSegmentIdBySegmentName(inputModel.GetName())
    tumorID = segmentationNode.GetSegmentation().GetSegmentIdBySegmentName("solidCavity")

    return segmentEditorWidget, segmentEditorNode, segmentationNode, breastID, tumorID


  def process(self, inputModel, probeTip, outputPoints, showResult=True):
    """
    Run the processing algorithm.s
    Can be used without GUI widget.
    :param inputModel: US probe 
    :param outputPoints: scanned point cloud
    :param showResult: show output volume in slice viewers
    """
    #self.pointCloud = self.generatePointCloud(outputPoints)
    self.model = self.generateModel(outputPoints, inputModel)

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

    

    """shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    exportFolderItemId = shNode.CreateFolderItem(shNode.GetSceneItemID(), "Cavity")

    segmentationNode.GetDisplayNode().SetSegmentVisibility(segmentID, 0)

    slicer.modules.segmentations.logic().ExportSegmentsToModels(segmentationNode, segmentID, exportFolderItemId)
    
    slicer.mrmlScene.RemoveNode(segmentationNode)"""




    


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
    logic.process(inputModel, outputPoints, threshold, True)
    outputScalarRange = outputPoints.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], threshold)

    # Test algorithm with inverted threshold
    logic.process(inputModel, outputPoints, threshold, False)
    outputScalarRange = outputPoints.GetImageData().GetScalarRange()
    self.assertEqual(outputScalarRange[0], inputScalarRange[0])
    self.assertEqual(outputScalarRange[1], inputScalarRange[1])

    self.delayDisplay('Test passed')
