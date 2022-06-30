"""
Python pour extraire une suite de jpeg pour ensuite ne faire un film avec makemovie.sh (ou faitfilm.sh)...
"""



########################################################################
# (appele dans fait_film_total.sh)
# Emploi : visit -s -cli -nowin vs_makemovie.py chemin_lata nom_champ chemin_movie nom_fichier_movie [first_state-laste_state]
#                                                   1          2              3                4            5
########################################################################

"""
chemin_lata = "/volatile/SIMU_STAGE/CUBE/SOURCE_QDM_GR/81_correction_QdM_finale_moyenne_glissante/terme_force_init_non_nul/avec_u_liq/signes_ok/RUN07/LATAS/cube.lata"
chemin_movie = "/volatile/SIMU_STAGE/CUBE/SOURCE_QDM_GR/81_correction_QdM_finale_moyenne_glissante/terme_force_init_non_nul/avec_u_liq/signes_ok/RUN07/LATAS"
nom_champ_lata = "VELOCITY_X_ELEM_DOM"
nom_champ = nom_champ_lata.rstrip("_DOM").rstrip("_ELEM").lower()
chemin_movie = "./MOVIE/"
nom_fichier_movie = "movie"
first_state,laste_state = 2,9
Lx,Ly,Lz =  .00259999999999999974 , .00259999999999999974 , .00259999999999999974 

"""

import sys
import os
import DNSTools3 as dtool

print(sys.argv)

# User gives the data, Reading the data
n_arg = 1
chemin_lata = sys.argv[n_arg]; n_arg+=1
print("chemin_lata : ",chemin_lata)
nom_champ_lata = sys.argv[n_arg]; n_arg+=1
nom_champ = nom_champ_lata.rstrip("_DOM").rstrip("_ELEM").lower()
print("nom_champ : ",nom_champ)
chemin_movie = sys.argv[n_arg]; n_arg+=1
print("chemin_movie : ",chemin_movie)
nom_fichier_movie = sys.argv[n_arg]; n_arg+=1
print("nom_fichier_movie : ",nom_fichier_movie)
first_state,laste_state = int(sys.argv[n_arg].lstrip("[").rstrip("]").split(",")[0]),int(sys.argv[n_arg].lstrip("[").rstrip("]").split(",")[-1])
print(first_state,laste_state)

# Domaine size
jdd = chemin_lata.split("/")[-1].replace(".lata",".data")
Lx = dtool.getParam(jdd, "uniform_domain_size_i")
Ly = dtool.getParam(jdd, "uniform_domain_size_j")
Lz = dtool.getParam(jdd, "uniform_domain_size_k")

# SetWindowLayout(2)
# SetActiveWindow(1)

# TimeSlider
curseur = CreateAnnotationObject("TimeSlider")
curseur.active = 1
curseur.position = (0.01,0.01)
curseur.width = 0.20
curseur.height = 0.035
curseur.text = 't = $time s'
curseur.timeFormatString = "%9.2f"

# Expresions
DefineScalarExpression(nom_champ, nom_champ_lata)
DefineScalarExpression("interfaces", "INTERFACES")

# Annotations
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.timeInfoFlag = 0
AnnotationAtts.legendInfoFlag = 1
AnnotationAtts.backgroundColor = (255, 255, 255, 255)
AnnotationAtts.foregroundColor = (0, 0, 0, 255)
AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.bboxFlag = 0
if not(AnnotationAtts.axes3D.visible == 0):
    AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Outside  # Inside, Outside, Both
    AnnotationAtts.axes3D.xAxis.tickMarks.visible = 0
    AnnotationAtts.axes3D.yAxis.tickMarks.visible = 0
    AnnotationAtts.axes3D.zAxis.tickMarks.visible = 0
    AnnotationAtts.axes3D.xAxis.title.visible = 0
    AnnotationAtts.axes3D.yAxis.title.visible = 0
    AnnotationAtts.axes3D.zAxis.title.visible = 0
    AnnotationAtts.axes3D.xAxis.title.userTitle = 1
    AnnotationAtts.axes3D.yAxis.title.userTitle = 1
    AnnotationAtts.axes3D.zAxis.title.userTitle = 1
    AnnotationAtts.axes3D.xAxis.title.title = "X"
    AnnotationAtts.axes3D.yAxis.title.title = "Y"
    AnnotationAtts.axes3D.zAxis.title.title = "Z"
    AnnotationAtts.axes3D.xAxis.title.font.scale = 1
    AnnotationAtts.axes3D.yAxis.title.font.scale = 1
    AnnotationAtts.axes3D.zAxis.title.font.scale = 1
    AnnotationAtts.axes3D.xAxis.label.font.scale = 2
    AnnotationAtts.axes3D.yAxis.label.font.scale = 2
    AnnotationAtts.axes3D.zAxis.label.font.scale = 2
    AnnotationAtts.axes3D.xAxis.label.visible = 0
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.axes3D.zAxis.label.visible = 0

SetAnnotationAttributes(AnnotationAtts)

# Contour Attributes
ContourAtts = ContourAttributes()
ContourAtts.defaultPalette.smoothing = ContourAtts.defaultPalette.None  # None, Linear, CubicSpline
ContourAtts.defaultPalette.equalSpacingFlag = 1
ContourAtts.defaultPalette.discreteFlag = 1
ContourAtts.defaultPalette.categoryName = "Standard"
ContourAtts.changedColors = ()
ContourAtts.colorType = ContourAtts.ColorByColorTable  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
ContourAtts.colorTableName = "difference"
ContourAtts.legendFlag = 1
ContourAtts.SetMultiColor(0, (255, 0, 0, 255))
ContourAtts.SetMultiColor(1, (0, 255, 0, 255))
ContourAtts.SetMultiColor(2, (0, 0, 255, 255))
ContourAtts.contourPercent = (10, 35, 65, 90)
ContourAtts.contourMethod = ContourAtts.Percent  # Level, Value, Percent
ContourAtts.minFlag = 1
ContourAtts.maxFlag = 1
if "lambda" in nom_champ:
    ContourAtts.min = -110000
    ContourAtts.max =  110000
else:
    ContourAtts.min = -0.1
    ContourAtts.max = 0.3
ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
ContourAtts.wireframe = 0

# Mesh Attributes - Interfaces
MeshAtts = MeshAttributes()
MeshAtts.legendFlag = 0
MeshAtts.lineWidth = 0
MeshAtts.meshColor = (180, 180, 180, 255)
MeshAtts.meshColorSource = MeshAtts.Foreground  # Foreground, MeshCustom, MeshRandom
MeshAtts.opaqueColorSource = MeshAtts.OpaqueCustom  # Background, OpaqueCustom, OpaqueRandom
MeshAtts.opaqueMode = MeshAtts.On  # Auto, On, Off
MeshAtts.pointSize = 0.05
MeshAtts.opaqueColor = (153, 153, 153, 255)
MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
MeshAtts.pointSizeVarEnabled = 0
MeshAtts.pointSizeVar = "default"
MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
MeshAtts.showInternal = 0
MeshAtts.pointSizePixels = 2
MeshAtts.opacity = 0.607843

# Large Box
LargeBoxAtts = BoxAttributes()
LargeBoxAtts.amount = LargeBoxAtts.Some  # Some, All
LargeBoxAtts.minx = -Lx
LargeBoxAtts.maxx = Lx # 0.007925 # 0.0317/4
LargeBoxAtts.miny = -Ly
LargeBoxAtts.maxy = Ly # 0.007925 # 0.0317/4
LargeBoxAtts.minz = -Lz
LargeBoxAtts.maxz = Lz # 0.007925 # 0.0317/4
LargeBoxAtts.inverse = 0
# Smoky Large Box
SubsetAtts_SLDOM = SubsetAttributes()
SubsetAtts_SLDOM.colorType = SubsetAtts_SLDOM.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
SubsetAtts_SLDOM.legendFlag = 0
SubsetAtts_SLDOM.singleColor = (0, 0, 0, 130)
SubsetAtts_SLDOM.subsetNames = ("Whole mesh (DOM)")
SubsetAtts_SLDOM.singleColor = (200, 200, 200, 255)
SubsetAtts_SLDOM.wireframe = 1
SubsetAtts_SLDOM.lineStyle = SubsetAtts_SLDOM.SOLID  # SOLID, DASH, DOT, DOTDASH
SubsetAtts_SLDOM.lineWidth = 2
SubsetAtts_SLDOM.opacity = 0.999

# Small Box
SmallBoxAtts = BoxAttributes()
SmallBoxAtts.amount = SmallBoxAtts.Some  # Some, All
SmallBoxAtts.minx = 0        
SmallBoxAtts.maxx = Lx/4 # 0.007925 # 0.0317/4
SmallBoxAtts.miny = 0       
SmallBoxAtts.maxy = Ly/4 # 0.007925 # 0.0317/4
SmallBoxAtts.minz = 0       
SmallBoxAtts.maxz = Lz/4 # 0.007925 # 0.0317/4
SmallBoxAtts.inverse = 0
# Showing Small Box
SubsetAtts_RDOM = SubsetAttributes()
SubsetAtts_RDOM.subsetNames = ("Zoom box")
SubsetAtts_RDOM.colorType = SubsetAtts_RDOM.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
SubsetAtts_RDOM.filledFlag = 1
SubsetAtts_RDOM.legendFlag = 0
SubsetAtts_RDOM.singleColor = (10, 180, 0, 255)
SubsetAtts_RDOM.opacity = 1
SubsetAtts_RDOM.wireframe = 1
SubsetAtts_RDOM.lineStyle = SubsetAtts_RDOM.SOLID  # SOLID, DASH, DOT, DOTDASH
SubsetAtts_RDOM.lineWidth = 3
# Translate and widen the small box
TransformAtts = TransformAttributes()
TransformAtts.doScale = 1
TransformAtts.scaleX = 3
TransformAtts.scaleY = 3
TransformAtts.scaleZ = 3
TransformAtts.doTranslate = 1
TransformAtts.translateX = Lx*1.02
TransformAtts.translateY = -(SmallBoxAtts.maxy / 2 )*TransformAtts.scaleX



# Opening latafile
OpenDatabase(chemin_lata, 0)

# Adjusting time
if first_state<0:
    first_state = 0
if laste_state<0:
    laste_state = TimeSliderGetNStates()
print("first_state,laste_state : ",first_state,laste_state)

# Plotting data
# - Contour
AddPlot("Contour", nom_champ, 1, 0)
# AddOperator("Box", 0)
# SetOperatorOptions(LargeBoxAtts, 0)
SetPlotOptions(ContourAtts)
# - Bubbles
AddPlot("Mesh", "interfaces", 1, 0)
# AddOperator("Box", 0)
# SetOperatorOptions(LargeBoxAtts, 0)
SetPlotOptions(MeshAtts)
# - The Small Box
AddPlot("Subset", "DOM", 0, 0)
AddOperator("Box", 0)
SetOperatorOptions(SmallBoxAtts, 0)
SetPlotOptions(SubsetAtts_RDOM)
# - The Smoky Large Box
AddPlot("Subset", "DOM", 0, 0)
SetPlotOptions(SubsetAtts_SLDOM)

# SetActiveWindow(2)
# SetActivePlots(0)
# DeleteActivePlots()
# SetActivePlots(0)
# DeleteActivePlots()
# SetActivePlots(0)
# DeleteActivePlots()

# Plotting data
# - Contour
AddPlot("Contour", nom_champ, 1, 0)
ContourAtts.legendFlag = 0
SetPlotOptions(ContourAtts)
AddOperator("Box", 0)
AddOperator("Transform", 0)
SetOperatorOptions(SmallBoxAtts, 2)
SetOperatorOptions(TransformAtts, 0)
# - Bubbles
AddPlot("Mesh", "interfaces", 1, 0)
SetPlotOptions(MeshAtts)
AddOperator("Box", 0)
AddOperator("Transform", 0)
SetOperatorOptions(SmallBoxAtts, 2)
SetOperatorOptions(TransformAtts, 0)
# - The Small Box
AddPlot("Subset", "DOM", 1, 0)
AddOperator("Box", 0)
AddOperator("Transform", 0)
SetOperatorOptions(SmallBoxAtts, 2)
SetOperatorOptions(TransformAtts, 0)
SetPlotOptions(SubsetAtts_RDOM)

#####################################
DrawPlots()
#####################################

# Hide MinMax
nom_champ_0 =  GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(nom_champ_0)
legend.drawMinMax = 0

ResetView()

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0, 0, 1)
View3DAtts.focus = (0.001665, -7.5e-05, -0.000897063)
View3DAtts.viewUp = (0, 1, 0)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 0.00744734
View3DAtts.nearPlane = -0.0148947
View3DAtts.farPlane = 0.0148947
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1.8
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (0.002665, -0.000975, -0.000897063)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
# End spontaneous state
# Fix Spontaneous state
ToggleMaintainViewMode()
 
# SetActiveWindow(1)
for i in range(first_state,laste_state):
    SetTimeSliderState(i)
    
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputDirectory = chemin_movie
    SaveWindowAtts.fileName = chemin_movie+nom_fichier_movie
    SaveWindowAtts.format = SaveWindowAtts.JPEG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
    SaveWindowAtts.width = 1050
    SaveWindowAtts.height = 980
    SaveWindowAtts.quality = 80
    SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate, LZW
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.pixelData = 1
    SaveWindowAtts.opts.help = ""
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
    
# SetActiveWindow(2)
# for i in range(first_state,laste_state):
    # SetTimeSliderState(i)
    # SaveWindowAtts = SaveWindowAttributes()
    # SaveWindowAtts.outputDirectory = chemin_movie
    # SaveWindowAtts.fileName = chemin_movie+"small_"+nom_fichier_movie
    # SaveWindowAtts.format = SaveWindowAtts.JPEG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
    # SaveWindowAtts.width = 1050
    # SaveWindowAtts.height = 980
    # SaveWindowAtts.quality = 80
    # SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate, LZW
    # SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
    # SaveWindowAtts.pixelData = 1
    # SaveWindowAtts.opts.help = ""
    # SetSaveWindowAttributes(SaveWindowAtts)
    # SaveWindow()
exit()
