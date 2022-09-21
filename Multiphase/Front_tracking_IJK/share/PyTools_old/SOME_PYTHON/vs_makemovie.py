"""
Python pour extraire une suite de jpeg pour ensuite ne faire un film avec makemovie.sh (ou faitfilm.sh)...
"""

import sys
import os

########################################################################
# (appele dans fait_film_total.sh)
# Emploi : visit -s -cli -nowin vs_makemovie.py chemin_lata nom_champ chemin_movie nom_fichier_movie [first_state-laste_state]
#                                                   1          2              3                4            5
########################################################################

"""
chemin_lata = "/volatile/SIMU_STAGE/CUBE/SOURCE_QDM_GR/81_correction_QdM_finale_moyenne_glissante/terme_force_init_non_nul/avec_u_liq/signes_ok/RUN07/LATAS/cube.lata"
chemin_movie = "/volatile/SIMU_STAGE/CUBE/SOURCE_QDM_GR/81_correction_QdM_finale_moyenne_glissante/terme_force_init_non_nul/avec_u_liq/signes_ok/RUN07/LATAS"
nom_champ = "VELOCITY_X_ELEM_DOM"
nom_fichier_movie = "movie0000"
first_state,laste_state = 2,9
"""

print(sys.argv)
# User gives the data
n_arg = 1
chemin_lata = sys.argv[n_arg]; n_arg+=1
print("chemin_lata : ",chemin_lata)
nom_champ = sys.argv[n_arg]; n_arg+=1
print("nom_champ : ",nom_champ)
chemin_movie = sys.argv[n_arg]; n_arg+=1
print("chemin_movie : ",chemin_movie)
nom_fichier_movie = sys.argv[n_arg]; n_arg+=1
print("nom_fichier_movie : ",nom_fichier_movie)
first_state,laste_state = int(sys.argv[n_arg].lstrip("[").rstrip("]").split("-")[0]),int(sys.argv[n_arg].lstrip("[").rstrip("]").split("-")[-1])
print("first_state,laste_state : ",first_state,laste_state)

# Opening latafile
OpenDatabase(chemin_lata, 0)

# Contour
AddPlot("Contour", nom_champ, 1, 1)
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
ContourAtts.min = -110000
ContourAtts.max = 110000
ContourAtts.scaling = ContourAtts.Linear  # Linear, Log
ContourAtts.wireframe = 1
SetPlotOptions(ContourAtts)

# Interfaces
AddPlot("Mesh", "INTERFACES", 1, 1)
MeshAtts = MeshAttributes()
MeshAtts.legendFlag = 1
MeshAtts.lineWidth = 0
MeshAtts.meshColor = (0, 0, 0, 255)
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
SetPlotOptions(MeshAtts)

DrawPlots()
SetActiveWindow(1)
#Source("/volatile/TRUST/trust-code/exec/VisIt/3.1.1/linux-x86_64/bin/makemovie.py")
ToggleCameraViewMode()

for i in range(first_state,laste_state):
    SetTimeSliderState(i+1)
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

exit()
