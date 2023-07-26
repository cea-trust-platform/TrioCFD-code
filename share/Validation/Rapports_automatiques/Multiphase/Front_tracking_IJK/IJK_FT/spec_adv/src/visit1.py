####################################
# Python pour obtenir des lineout avec visit
####################################


from sys import *
import os

# CHEMIN
pwd = os.getcwd()
print("pwd",pwd)

# NOM DU CAS
cas = pwd.split("/")[-2]
print("cas",cas)

# EMPLACEMENT DE LA FICHE (build)
fiche = pwd.replace("/"+cas+"/RUN00","")

L = 0.004
time_iterations = [1,6,11,16,21]
lineout = [(-L/2., 0, 0),(L/2, 0, 0)]

# PARAMETRAGE EN FONCTION DU CAS
if ("X00" in cas) or ("100" in cas):
    which_direction = "X" 
    time_iterations = [1,10,19,26]
    only_u = True
elif ("0X0" in cas) or ("010" in cas):
    which_direction = "Y" 
    if ("0X0" in cas):
        time_iterations = [1,6,11,16,17]
    if ("010" in cas):
        lineout = [(0, -L/2., 0),(0, L/2, 0)]
        time_iterations = [1,6,12,18]
elif ("001" in cas):
    which_direction = "Z"
    lineout = [(0, 0, -L/2.),(0, 0, L/2)]
    time_iterations = [1,6,12,18]
elif ("101" in cas):
    exit()
print("which_direction",which_direction)
################################################################

##### Commandes Visit ###############################################
# OpenDatabase(pwd+"/build/"+fiche+"/RUN00/LATAS/spec_bulles_point2.lata", 0)
OpenDatabase(pwd+"/LATAS/spec_bulles_point2.lata", 0)
# OpenDatabase("localhost:/volatile/FFTW/THI_FFT/Rapport_validation/spec_adv/build/ADV_X00/RUN00/LATAS/spec_bulles_point2.lata", 0)

AddPlot("Curve", "operators/Lineout/PRESSURE_ELEM_DOM", 1, 1)
AddPlot("Curve", "operators/Lineout/FORCE_PH_"+which_direction+"_FACES_DOM_dual", 1, 1)
AddPlot("Curve", "operators/Lineout/VELOCITY_"+which_direction+"_FACES_DOM_dual", 1, 1)

SetActivePlots((0, 1))
SetActivePlots(0)
CurveAtts = CurveAttributes()
CurveAtts.showLines = 1
CurveAtts.lineWidth = 3
CurveAtts.curveColorSource = CurveAtts.Custom  # Cycle, Custom
CurveAtts.curveColor = (50, 255, 50, 255)
CurveAtts.showLegend = 1
CurveAtts.showLabels = 0
SetPlotOptions(CurveAtts)

SetActivePlots((0, 1))
SetActivePlots(1)
CurveAtts = CurveAttributes()
CurveAtts.showLines = 1
CurveAtts.lineWidth = 3
CurveAtts.showLegend = 1
CurveAtts.showLabels = 0
CurveAtts.curveColorSource = CurveAtts.Custom  # Cycle, Custom
CurveAtts.curveColor = (255, 0, 0, 255)
SetPlotOptions(CurveAtts)

SetActivePlots((0, 1))
SetActivePlots(2)
CurveAtts = CurveAttributes()
CurveAtts.showLines = 1
CurveAtts.lineWidth = 3
CurveAtts.showLegend = 1
CurveAtts.showLabels = 0
CurveAtts.curveColorSource = CurveAtts.Custom  # Cycle, Custom
CurveAtts.curveColor = (50, 50, 255, 255)
SetPlotOptions(CurveAtts)

LineoutAtts = LineoutAttributes()
LineoutAtts.point1 = lineout[0] #(-L/2., 0, 0)
LineoutAtts.point2 = lineout[1] #(L/2, 0, 0)
LineoutAtts.numberOfSamplePoints = 50
SetOperatorOptions(LineoutAtts, 0, 1)
SetActivePlots((0, 1))

# Begin spontaneous state
ViewCurveAtts = ViewCurveAttributes()
ViewCurveAtts.domainCoords = (0, L)
ViewCurveAtts.rangeCoords = (-2,2)
ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
ViewCurveAtts.domainScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
ViewCurveAtts.rangeScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
SetViewCurve(ViewCurveAtts)
# End spontaneous state

for t in time_iterations:
    SetTimeSliderState(t)
    DrawPlots()
    
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 0
    SaveWindowAtts.outputDirectory = fiche 
    SaveWindowAtts.fileName = "lineout_"+cas+"_t"+str(t)
    SaveWindowAtts.family = 0
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
    SaveWindowAtts.width = 1024
    SaveWindowAtts.height = 1024
    SaveWindowAtts.quality = 80
    # #SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate, LZW
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.pixelData = 1
    SaveWindowAtts.subWindowAtts.win1.position = (0, 0)
    SaveWindowAtts.subWindowAtts.win1.size = (128, 128)
    SaveWindowAtts.opts.types = ()
    SaveWindowAtts.opts.help = ""
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()

exit()
