#!/usr/bin/python

import os, shutil, glob
from math import pi

l = glob.glob('DNS.lata.*.coordx')
l.extend(glob.glob('latas/DNS.lata.*.coordx'))
dom=l[0].split('/')[-1].replace('DNS.lata.', "").replace(".coordx", "").replace("_EXT", "")
#dom="DOM"
#if (os.path.isfile('latas/DNS.lata.DOM_NEW_EXT.coordx') or os.path.isfile('DNS.lata.DOM_NEW_EXT.coordx')):
#   dom="DOM_NEW"
#elif (os.path.isfile('latas/DNS.lata.DOM_EXT.coordx') or os.path.isfile('DNS.lata.DOM_EXT.coordx')):
#   dom="DOM"
#else:
#   raise Exception("DOM or DOM_NEW mesh not found")
#   pass

def defineLambda2():
  DefineScalarExpression("lambda2", "LAMBDA2_ELEM_%s"%dom)
  return

def defineLambda2Visit():
  DefineVectorExpression("velocity", "{VELOCITY_X_FACES_%s_dual,VELOCITY_Y_FACES_%s_dual,VELOCITY_Z_FACES_%s_dual}"%(dom,dom,dom))
  DefineScalarExpression("vorticity", "magnitude(curl(velocity))")
  DefineScalarExpression("du", "gradient(VELOCITY_X_FACES_%s_dual)"%dom)
  DefineScalarExpression("dv", "gradient(VELOCITY_Y_FACES_%s_dual)"%dom)
  DefineScalarExpression("dw", "gradient(VELOCITY_Z_FACES_%s_dual)"%dom)
  DefineScalarExpression("L11", "du[2]*du[2] + du[1]*dv[2] + du[0]*dw[2]")
  DefineScalarExpression("L12", "((du[1] + dv[2])* (du[2] + dv[1]) + dv[0]*dw[2] + du[0]*dw[1])/2.0")
  DefineScalarExpression("L13", "(dv[2]*dv[0] + du[1]*dw[1] + (du[0]+dw[2])*(du[2]+dw[0]))/2.0")
  DefineScalarExpression("L22", "du[1]*dv[2] + dv[1]* dv[1] + dv[0]*dw[1]")
  DefineScalarExpression("L23", "(du[0]*dv[2]+du[1]*dw[2] + (dv[0]+dw[1])*(dv[1]+ dw[0]))/2.0")
  DefineScalarExpression("L33", "du[0]*dw[2] + dv[0]*dw[1]+ dw[0]*dw[0]")
  DefineScalarExpression("lambda", "eigenvalue({{L11,L12,L13}, {L12,L22,L23}, {L13,L23,L33}})")
  DefineScalarExpression("lambda2", " lambda[1]")
  return


DefineScalarExpression("ln_curvature", "ln(max(0.01, -COURBURE_SOM_INTERFACES))")
DefineScalarExpression("curvature", "COURBURE_SOM_INTERFACES") # So that it has a nice name shown!
DefineScalarExpression("temperature", "TEMPERATURE_0_ELEM_%s"%dom) # So that it has a nice name shown!
opts="GB" 
opts="ABC" # New colors version from Antoine

# Defining some attributes : 
if (opts == "ABC" ) : 
     ContourAtts = ContourAttributes()
     ContourAtts.defaultPalette.GetControlPoints(0).colors = (255, 0, 0,255)
     ContourAtts.defaultPalette.GetControlPoints(0).position = 0
     ContourAtts.defaultPalette.GetControlPoints(1).colors = (0, 255, 0,255)
     ContourAtts.defaultPalette.GetControlPoints(1).position = 0.034
     ContourAtts.defaultPalette.GetControlPoints(2).colors = (0, 0, 255,255)
     ContourAtts.defaultPalette.GetControlPoints(2).position = 0.069
     ContourAtts.defaultPalette.GetControlPoints(3).colors = (0, 255, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(3).position = 0.103
     ContourAtts.defaultPalette.GetControlPoints(4).colors = (255, 0, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(4).position = 0.138
     ContourAtts.defaultPalette.GetControlPoints(5).colors = (255, 255, 0, 255)
     ContourAtts.defaultPalette.GetControlPoints(5).position = 0.172
     ContourAtts.defaultPalette.GetControlPoints(6).colors = (255, 135, 0, 255)
     ContourAtts.defaultPalette.GetControlPoints(6).position = 0.207
     ContourAtts.defaultPalette.GetControlPoints(7).colors = (255, 0, 135, 255)
     ContourAtts.defaultPalette.GetControlPoints(7).position = 0.241
     ContourAtts.defaultPalette.GetControlPoints(8).colors = (168, 168, 168, 255)
     ContourAtts.defaultPalette.GetControlPoints(8).position = 0.276
     ContourAtts.defaultPalette.GetControlPoints(9).colors = (255, 68, 68, 255)
     ContourAtts.defaultPalette.GetControlPoints(9).position = 0.31
     ContourAtts.defaultPalette.GetControlPoints(10).colors = (99, 255, 99, 255)
     ContourAtts.defaultPalette.GetControlPoints(10).position = 0.345
     ContourAtts.defaultPalette.GetControlPoints(11).colors = (99, 99, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(11).position = 0.379
     ContourAtts.defaultPalette.GetControlPoints(12).colors = (40, 165, 165, 255)
     ContourAtts.defaultPalette.GetControlPoints(12).position = 0.414
     ContourAtts.defaultPalette.GetControlPoints(13).colors = (255, 99, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(13).position = 0.448
     ContourAtts.defaultPalette.GetControlPoints(14).colors = (255, 255, 99, 255)
     ContourAtts.defaultPalette.GetControlPoints(14).position = 0.483
     ContourAtts.defaultPalette.GetControlPoints(15).colors = (255, 170, 99, 255)
     ContourAtts.defaultPalette.GetControlPoints(15).position = 0.517
     ContourAtts.defaultPalette.GetControlPoints(16).colors = (170, 79, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(16).position = 0.552
     ContourAtts.defaultPalette.GetControlPoints(17).colors = (150, 0, 0, 255)
     ContourAtts.defaultPalette.GetControlPoints(17).position = 0.586
     ContourAtts.defaultPalette.GetControlPoints(18).colors = (0, 150, 0, 255)
     ContourAtts.defaultPalette.GetControlPoints(18).position = 0.621
     ContourAtts.defaultPalette.GetControlPoints(19).colors = (0, 0, 150, 255)
     ContourAtts.defaultPalette.GetControlPoints(19).position = 0.655
     ContourAtts.defaultPalette.GetControlPoints(20).colors = (0, 109, 109, 255)
     ContourAtts.defaultPalette.GetControlPoints(20).position = 0.69
     ContourAtts.defaultPalette.GetControlPoints(21).colors = (150, 0, 150, 255)
     ContourAtts.defaultPalette.GetControlPoints(21).position = 0.724
     ContourAtts.defaultPalette.GetControlPoints(22).colors = (150, 150, 0, 255)
     ContourAtts.defaultPalette.GetControlPoints(22).position = 0.759
     ContourAtts.defaultPalette.GetControlPoints(23).colors = (150, 84, 0, 255)
     ContourAtts.defaultPalette.GetControlPoints(23).position = 0.793
     ContourAtts.defaultPalette.GetControlPoints(24).colors = (160, 0, 79, 255)
     ContourAtts.defaultPalette.GetControlPoints(24).position = 0.828
     ContourAtts.defaultPalette.GetControlPoints(25).colors = (255, 104, 28, 255)
     ContourAtts.defaultPalette.GetControlPoints(25).position = 0.862
     ContourAtts.defaultPalette.GetControlPoints(26).colors = (0, 170, 81, 255)
     ContourAtts.defaultPalette.GetControlPoints(26).position = 0.897
     ContourAtts.defaultPalette.GetControlPoints(27).colors = (68, 255, 124, 255)
     ContourAtts.defaultPalette.GetControlPoints(27).position = 0.931
     ContourAtts.defaultPalette.GetControlPoints(28).colors = (0, 130, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(28).position = 0.966
     ContourAtts.defaultPalette.GetControlPoints(29).colors = (130, 0, 255, 255)
     ContourAtts.defaultPalette.GetControlPoints(29).position = 1
     ContourAtts.defaultPalette.smoothing = ContourAtts.defaultPalette.None  # None, Linear, CubicSpline
     ContourAtts.defaultPalette.equalSpacingFlag = 1
     ContourAtts.defaultPalette.discreteFlag = 1
     ContourAtts.defaultPalette.categoryName = "Standard"
     ContourAtts.changedColors = (1, 0, 2)
     ContourAtts.colorType = ContourAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
     ContourAtts.colorTableName = "Default"
     ContourAtts.invertColorTable = 0
     ContourAtts.legendFlag = 1
     ContourAtts.lineStyle = ContourAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
     ContourAtts.lineWidth = 9
     ContourAtts.singleColor = (0, 128, 128, 52)
     ContourAtts.SetMultiColor(0, (255, 153, 0, 193))
     ContourAtts.SetMultiColor(1, (153, 204, 0, 129))
     ContourAtts.SetMultiColor(2, (51, 204, 204, 66))
     ContourAtts.SetMultiColor(3, (0, 255, 255, 255))
     ContourAtts.SetMultiColor(4, (255, 0, 255, 255))
     ContourAtts.SetMultiColor(5, (255, 255, 0, 255))
     ContourAtts.SetMultiColor(6, (255, 135, 0, 255))
     ContourAtts.SetMultiColor(7, (255, 0, 135, 255))
     ContourAtts.SetMultiColor(8, (168, 168, 168, 255))
     ContourAtts.SetMultiColor(9, (255, 68, 68, 255))
     ContourAtts.contourNLevels = 3
     ContourAtts.contourValue = (-100, -10, -1)
     ContourAtts.contourPercent = ()
     ContourAtts.contourMethod = ContourAtts.Value  # Level, Value, Percent
     ContourAtts.minFlag = 1
     ContourAtts.maxFlag = 1
     ContourAtts.min = -100000
     ContourAtts.max = 0
     ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
     ContourAtts.wireframe = 0
else:
     ContourAtts = ContourAttributes()
     ContourAtts.contourValue = (-1)
     ContourAtts.contourPercent = ()
     ContourAtts.contourMethod = ContourAtts.Value  # Level, Value, Percent
     # Palette : 
     ContourAtts.changedColors = (0)
     ContourAtts.colorType = ContourAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
     ContourAtts.colorTableName = "Default"
     ContourAtts.invertColorTable = 0
     ContourAtts.legendFlag = 0 # Masquer cette legende
     ContourAtts.lineStyle = ContourAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
     ContourAtts.lineWidth = 0
     ContourAtts.singleColor = (128, 128, 128, 110)
     pass


dest="Figs"
frequency = 1
opt=1 # With lambda2
# opt=0 # Without Lambda2
showL2inChannel=False
#showL2inChannel=True
optZoom = 1 # With the zoomed box
# optZoom = 0 # Without the zoomed box.
visu = "3D"
visu = "2D"

if (opts=="ABC"):
   kappa_min = -40
   kappa_max = 2
else:
   kappa_min = 0
   kappa_max = 3
   pass

if (os.path.isfile("./latas/compil.lata")):
   inp="./latas/compil.lata"
elif (os.path.isfile("DNS.lata")):
   inp="DNS.lata"
else:
   raise Exception("input lata File DNS.lata or ./latas/compil.lata is missing")

if (visu == "3D"): dest+="3D"
if opt: dest+=opts+"_Lambda2"

OpenDatabase(inp, 0)
if opt: defineLambda2()

with_bubbles=True
realBubblesOnly=1
try:
   if (opts=="ABC"):
      AddPlot("Pseudocolor", "curvature", 0, 0)
      vname="curvature"
   else:
      AddPlot("Pseudocolor", "ln_curvature", 0, 0)
      vname="curvature (ln)"
      pass
except:
   # The curvature is not present, we assume it's a single-phase calculation. 
   with_bubbles=False
   pass

# Changing parameter for the curvature:
if with_bubbles:
   PseudocolorAtts = PseudocolorAttributes()
   PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
   PseudocolorAtts.skewFactor = 1
   PseudocolorAtts.limitsMode = PseudocolorAtts.CurrentPlot  # OriginalData, CurrentPlot
   PseudocolorAtts.minFlag = 1 # GB
   PseudocolorAtts.min = kappa_min
   PseudocolorAtts.maxFlag = 1 # GB 
   PseudocolorAtts.max = kappa_max
   PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
   PseudocolorAtts.invertColorTable = 0
   PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
   PseudocolorAtts.opacityVariable = ""
   PseudocolorAtts.opacity = 1
   PseudocolorAtts.opacityVarMin = 0
   PseudocolorAtts.opacityVarMax = 1
   PseudocolorAtts.opacityVarMinFlag = 0
   PseudocolorAtts.opacityVarMaxFlag = 0
   PseudocolorAtts.pointSize = 0.05
   PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
   PseudocolorAtts.pointSizeVarEnabled = 0
   PseudocolorAtts.pointSizeVar = "default"
   PseudocolorAtts.pointSizePixels = 2
   PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
   PseudocolorAtts.lineWidth = 0
   PseudocolorAtts.tubeResolution = 10
   PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
   PseudocolorAtts.tubeRadiusAbsolute = 0.125
   PseudocolorAtts.tubeRadiusBBox = 0.005
   PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox # Absolute, FractionOfBBox
   PseudocolorAtts.renderSurfaces = 1
   PseudocolorAtts.renderWireframe = 0
   PseudocolorAtts.renderPoints = 0
   PseudocolorAtts.smoothingLevel = 0
   PseudocolorAtts.legendFlag = 1
   PseudocolorAtts.lightingFlag = 1
   if (opts=="ABC"):
      PseudocolorAtts.colorTableName = "gray"
      PseudocolorAtts.tubeRadiusVarEnabled = 0
      PseudocolorAtts.tubeRadiusVar = ""
      PseudocolorAtts.tubeRadiusVarRatio = 10
      PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
      PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
      PseudocolorAtts.endPointRadiusAbsolute = 0.125
      PseudocolorAtts.endPointRadiusBBox = 0.05
      PseudocolorAtts.endPointResolution = 10
      PseudocolorAtts.endPointRatio = 5
      PseudocolorAtts.endPointRadiusVarEnabled = 0
      PseudocolorAtts.endPointRadiusVar = ""
      PseudocolorAtts.endPointRadiusVarRatio = 10
   else:
      PseudocolorAtts.colorTableName = "bluehot" # "PuBu"
      PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
      PseudocolorAtts.endPointRadiusAbsolute = 1
      PseudocolorAtts.endPointRadiusBBox = 0.005
      PseudocolorAtts.endPointRatio = 2
      pass
   SetPlotOptions(PseudocolorAtts)
   realBubblesOnly=1
   if realBubblesOnly:
      AddOperator("Threshold", 0)
      ThresholdAtts = ThresholdAttributes()
      ThresholdAtts.outputMeshType = 0
      ThresholdAtts.listedVarNames = ("default", "COMPO_CONNEXE_ELEM_INTERFACES")
      ThresholdAtts.zonePortions = (1, 1)
      ThresholdAtts.lowerBounds = (-1e+37, -1e+37)
      ThresholdAtts.upperBounds = (1e+37, 1e+37)
      ThresholdAtts.defaultVarName = "default"
      ThresholdAtts.defaultVarIsScalar = 0
      SetOperatorOptions(ThresholdAtts, 0)
      print "Selecting real bubbles only for the curvature."
      pass
   DrawPlots()
   print "Ploting the curvature on the interfaces and changing rendering properties"
   pass

BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.Some  # Some, All
BoxAtts.minx = 1 # pi
BoxAtts.maxx = 2.8 # pi*2.
BoxAtts.miny = 0 # pi*2./4.
BoxAtts.maxy = 1 # pi
BoxAtts.minz = 0
BoxAtts.maxz = 1
BoxAtts.inverse = 0

TransformAtts = TransformAttributes()
TransformAtts.doRotate = 0
TransformAtts.rotateOrigin = (0, 0, 0)
TransformAtts.rotateAxis = (0, 0, 1)
TransformAtts.rotateAmount = 0
TransformAtts.rotateType =    TransformAtts.Deg  # Deg, Rad
TransformAtts.doScale = 1
TransformAtts.scaleOrigin = (0, 0, 0)
TransformAtts.scaleX = 3
TransformAtts.scaleY = 3
TransformAtts.scaleZ = 3
TransformAtts.doTranslate = 1
TransformAtts.translateX = -3
TransformAtts.translateY = 0
TransformAtts.translateZ = 2.5
TransformAtts.transformType =    TransformAtts.Similarity  # Similarity, Coordinate, Linear
TransformAtts.inputCoordSys =    TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
TransformAtts.outputCoordSys =    TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
TransformAtts.continuousPhi = 0
TransformAtts.m00 = 1
TransformAtts.m01 = 0
TransformAtts.m02 = 0
TransformAtts.m03 = 0
TransformAtts.m10 = 0
TransformAtts.m11 = 1
TransformAtts.m12 = 0
TransformAtts.m13 = 0
TransformAtts.m20 = 0
TransformAtts.m21 = 0
TransformAtts.m22 = 1
TransformAtts.m23 = 0
TransformAtts.m30 = 0
TransformAtts.m31 = 0
TransformAtts.m32 = 0
TransformAtts.m33 = 1
TransformAtts.invertLinearTransform = 0
TransformAtts.vectorTransformMethod =    TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
TransformAtts.transformVectors = 1

# Trace du lambda2 dans le zoom
if (opt and optZoom): 
   AddPlot("Contour", "lambda2", 0, 0)
   #ChangeActivePlotsVar("COMPO_CONNEXE_ELEM_INTERFACES")
   AddOperator("Box", 0)
   SetActivePlots(1)
   SetOperatorOptions(BoxAtts, 0, 0) # Last option with 0 stands for 'not applied to all plots'
   #
   AddOperator("Transform", 0)
   SetOperatorOptions(TransformAtts, 0, 0)
   #
   SetPlotOptions(ContourAtts)
   DrawPlots()
   print "Ploting lambda2 in the shifted box"
   pass

AddPlot("Subset", dom+"_EXT", 1, 0) # ABC met 1, 0), je mets 0,0) why? 
SubsetAtts = SubsetAttributes()
SubsetAtts.colorType = SubsetAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
SubsetAtts.invertColorTable = 0
SubsetAtts.filledFlag = 1
SubsetAtts.legendFlag = 0 # Masquer cette legende
SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
SubsetAtts.lineWidth = 0
if ((opts == "GB") or 1 ):
   SubsetAtts.colorTableName = "Default"
   SubsetAtts.singleColor = (0, 0, 0, 14)
else:
   SubsetAtts.singleColor = (0, 0, 0, 255)
   SubsetAtts.SetMultiColor(0, (255, 0, 0, 0))
   pass

SubsetAtts.subsetNames = ("Whole mesh")
#SubsetAtts.subsetType = SubsetAtts.Mesh  # Domain, Group, Material, EnumScalar, Mesh, Unknown
SubsetAtts.opacity = 1
SubsetAtts.wireframe = 0
SubsetAtts.drawInternal = 0
SubsetAtts.smoothingLevel = 0
SubsetAtts.pointSize = 0.05
SubsetAtts.pointType = SubsetAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
SubsetAtts.pointSizeVarEnabled = 0
SubsetAtts.pointSizeVar = "default"
SubsetAtts.pointSizePixels = 2
SetPlotOptions(SubsetAtts)
DrawPlots()
print "Showing the shaded box for the channel"

AddPlot("Subset", dom+"_EXT", 0, 0)
SetActivePlots(2+opt)
SubsetAtts = SubsetAttributes()
SubsetAtts.colorType = SubsetAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
SubsetAtts.filledFlag = 1
SubsetAtts.legendFlag = 0
SubsetAtts.singleColor = (0, 0, 0, 255)
SubsetAtts.opacity = 1
SubsetAtts.wireframe = 1
SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
SubsetAtts.lineWidth = 1
SetPlotOptions(SubsetAtts)
print "Showing the contours (wireframe) of the main channel"
DrawPlots()

if optZoom:
   if (opts=="ABC"):
      AddPlot("Pseudocolor", "curvature", 0, 0)
      vname="curvature"
   else:
      AddPlot("Pseudocolor", "ln_curvature", 0, 0)
      vname="curvature (ln)"
      pass
   SetActivePlots(3+opt)
   AddOperator("Box", 0)
   AddOperator("Transform", 0)
   SetOperatorOptions(BoxAtts, 3+opt, 0) # Last option with 0 stands for 'not applied to all plots'
   SetOperatorOptions(TransformAtts, 0, 0)
   PseudocolorAtts.legendFlag = 0
   SetPlotOptions(PseudocolorAtts)
   print "Showing the curvature on interfaces in the shifted small box"
   DrawPlots()
   #
   AddPlot("Subset", dom+"_EXT", 0, 0)
   SetActivePlots(4+opt)
   AddOperator("Box", 0)
   AddOperator("Transform", 0)
   SetOperatorOptions(BoxAtts, 0, 0) # Last option with 0 stands for 'not applied to all plots'
   SetOperatorOptions(TransformAtts, 0, 0)
   SubsetAtts_RDOM = SubsetAttributes()
   SubsetAtts_RDOM.colorType = SubsetAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
   SubsetAtts_RDOM.filledFlag = 1
   SubsetAtts_RDOM.legendFlag = 0
   SubsetAtts_RDOM.singleColor = (255, 0, 0, 255)
   SubsetAtts_RDOM.opacity = 1
   SubsetAtts_RDOM.wireframe = 1
   SubsetAtts_RDOM.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
   SubsetAtts_RDOM.lineWidth = 3
   SetPlotOptions(SubsetAtts_RDOM)
   print "Showing the red contours (wireframe) of the small shifted box"
   DrawPlots()
   #
   if (not showL2inChannel): 
      # Otherwise, we cannot see a thing!
      AddPlot("Subset", dom+"_EXT", 0, 0)
      SetActivePlots(5+opt)
      AddOperator("Box", 0)
      SetOperatorOptions(BoxAtts, 0, 0) # Last option with 0 stands for "not applied to all plots"
      # SubsetAtts_RDOM.singleColor = (0, 0, 0, 255) # -> if we like it in black... 
      SetPlotOptions(SubsetAtts_RDOM)
      print "Showing the red contours (wireframe) of the small box (before transformation)"
      DrawPlots()
      pass
   pass

# Walls : 
for zwall in [0., 2.]:
   AddPlot("Subset", dom+"_EXT", 0, 0)
   AddOperator("Slice", 0)
   # SetActivePlots(4)
   SliceAtts = SliceAttributes()
   SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
   SliceAtts.originPoint = (0, 0, zwall)
   SliceAtts.normal = (0, 0, 1)
   SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
   SliceAtts.upAxis = (0, 1, 0)
   SliceAtts.project2d = 0
   SliceAtts.interactive = 1
   SliceAtts.flip = 0
   SetOperatorOptions(SliceAtts, 0)
   SubsetAtts_wall = SubsetAttributes()
   SubsetAtts_wall.colorType = SubsetAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
   SubsetAtts_wall.filledFlag = 1
   SubsetAtts_wall.legendFlag = 0
   SubsetAtts_wall.singleColor = (0, 0, 0, 80)
   SubsetAtts_wall.opacity = 1
   SetPlotOptions(SubsetAtts_wall)
   print "Showing walls on the main channel"
   DrawPlots()
   pass


# Trace du lambda2 dans la mainbox 
if (opt and showL2inChannel): 
   print "WARNING : showing lambda2 on the whole channel will take some time to compute isosurface... be ready for a coffea..."
   AddPlot("Contour", "lambda2", 0, 0)
   if (optZoom):
     # Lambda2 legend is already shown by zoom... We dont want it twice
     ContourAtts.legendFlag = 0
     pass
   SetPlotOptions(ContourAtts)
   DrawPlots()
   print "Plotting lambda2 in the main channel"
   pass

# Annotations : 
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes3D.visible = 1
AnnotationAtts.axes3D.autoSetTicks = 1
AnnotationAtts.axes3D.autoSetScaling = 1
AnnotationAtts.axes3D.lineWidth = 0
AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
AnnotationAtts.axes3D.triadFlag = 0
AnnotationAtts.axes3D.xAxis.title.visible = 0
AnnotationAtts.axes3D.xAxis.label.visible = 1
AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
AnnotationAtts.axes3D.xAxis.label.font.scale = 0.7
AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
AnnotationAtts.axes3D.xAxis.label.font.bold = 0
AnnotationAtts.axes3D.xAxis.label.font.italic = 0
AnnotationAtts.axes3D.xAxis.label.scaling = 0
AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 1
AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 2.8
AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.05
AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.5
AnnotationAtts.axes3D.xAxis.grid = 0
AnnotationAtts.axes3D.yAxis.title.visible = 0
AnnotationAtts.axes3D.yAxis.label.visible = 0
AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
AnnotationAtts.axes3D.yAxis.grid = 0
AnnotationAtts.axes3D.zAxis.title.visible = 0
AnnotationAtts.axes3D.zAxis.label.visible = 0
AnnotationAtts.axes3D.zAxis.tickMarks.visible = 0
AnnotationAtts.axes3D.bboxFlag = 1
AnnotationAtts.axes3D.setBBoxLocation = 1
AnnotationAtts.axes3D.bboxLocation = (1, 2.8, 0, 1, 0, 1)
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.timeInfoFlag = 1
AnnotationAtts.legendInfoFlag = 1
SetAnnotationAttributes(AnnotationAtts)

# TimeSlider : 
curseur = CreateAnnotationObject("TimeSlider")
curseur.active = 1
curseur.position = (0.01,0.01)
# ABC curseur.position = (0.01,0.2)
curseur.width = 0.20
curseur.height = 0.035
curseur.text = 't = $time s'
curseur.timeFormatString = "%9.2f"
# Text
curseur2 = CreateAnnotationObject("Text2D")
curseur2.text = "G. Bois & A. du Cluzeau \nCEA/DEN 2018-2019"
# curseur2.text = "TrioIJK \nG. Bois & A. du Cluzeau \nCEA/DEN 2018-2019"
curseur2.visible = 1
curseur2.active = 1
curseur2.position = (0.01,0.055)
curseur2.height = 0.02
curseur2.textColor = (0, 0, 0, 180)
curseur2.useForegroundForTextColor = 1
#curseur2.fontFamily = Arial  # Arial, Courier, Times
curseur2.fontBold = 0
curseur2.fontItalic = 0
curseur2.fontShadow = 0

# Image :
# logo du haut : 
curseur3 = CreateAnnotationObject("Image")
curseur3.image = ("~/logo_cea.png")
curseur3.width = 40
curseur3.position = (0.01, 0.20)
curseur3.visible = 1

# logo du dessous : 
curseur4 = CreateAnnotationObject("Image")
curseur4.image = ("~/logo_TrioIJK.png")
curseur4.width = 60
curseur4.position = (0.01, 0.10)
curseur4.visible = 1

ResetView()
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0, 1, 0)
View3DAtts.viewUp = (1, 0, 0)
View3DAtts.perspective = 1
View3DAtts.focus = (pi, pi/2, 1.)
View3DAtts.centerOfRotation = (pi, pi/2, 1.)

# Begin spontaneous state
View3DAtts.focus = (3.4, 1.54396, 2.4)
View3DAtts.viewNormal = (0, 1, -0.06)
View3DAtts.viewAngle = 20 # Joue sur l'angle entre les 2 murs...
View3DAtts.parallelScale = 4.25 #   anciennement 4.16022 
View3DAtts.nearPlane = -8.32044
View3DAtts.farPlane = 8.32044
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2 #   anciennement 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (3.32327, 1.54396, 1)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state

ViewAxisArrayAtts = ViewAxisArrayAttributes()
ViewAxisArrayAtts.domainCoords = (0, 1)
ViewAxisArrayAtts.rangeCoords = (0, 1)
ViewAxisArrayAtts.viewportCoords = (0.15, 0.9, 0.1, 0.85)
SetViewAxisArray(ViewAxisArrayAtts)

if (os.path.isdir(dest)) : shutil.rmtree(dest)
os.mkdir(dest)
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = dest
SaveWindowAtts.fileName = "simple_view"
SaveWindowAtts.family = 0      # To add left and right prefix... but also stops the automatic numbering 
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
# Full HD : 
if (visu == "3D"):
   SaveWindowAtts.width = 3360
   SaveWindowAtts.height = 1896
else:
   SaveWindowAtts.width = 1920
   SaveWindowAtts.height = 1080
   pass

SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
if (visu == "3D"):
   SaveWindowAtts.stereo = 1 # Pour la 3D... 
else:
   SaveWindowAtts.stereo = 0
   pass

SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint
SaveWindowAtts.advancedMultiWindowSave = 0
SaveWindowAtts.format = SaveWindowAtts.JPEG
SaveWindowAtts.quality = 100
SetSaveWindowAttributes(SaveWindowAtts)

# Fixer la vue : 
ToggleLockViewMode()
#import pdb; pdb.set_trace()

# Cache le minmax de la courbure
nom_champ =  GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(nom_champ)
legend.drawMinMax = 0
if (opt):
   # Cache le minmax du lambda2
   nom_champl2 =  GetPlotList().GetPlots(1).plotName
   legendl2 = GetAnnotationObject(nom_champl2)
   legendl2.drawMinMax = 0
   pass

# for state in range(0, 1):
count = 0
for state in range(0,TimeSliderGetNStates(), frequency):
   TimeSliderSetState(state)
   DrawPlots()
   SaveWindowAtts.fileName = "simple_view_%05d"%count
   SetSaveWindowAttributes(SaveWindowAtts)
   SaveWindow()
   count+=1
   pass
   
print "The end of figures generation!"
SaveSession("./my_movie.session")
exit()
