a = AnnotationAttributes()
r = RenderingAttributes()
s = SliceAttributes()
b = BoxAttributes()
t = TransformAttributes()
p = PseudocolorAttributes()
v = GetView2D()

# Ouverture du fichier
#---------------------

OpenDatabase("Cas.lata")
SetTimeSliderState(TimeSliderGetNStates()-1)

# Parametres du rendu
#--------------------

r.antialiasing = 1
SetRenderingAttributes(r)

# Parametres des annotations
#---------------------------

a.databaseInfoFlag = 0
a.timeInfoFlag = 1
a.userInfoFlag = 0
a.legendInfoFlag = 1
a.axes2D.visible = 0
SetAnnotationAttributes(a)

# Plot 0
#-------

#AddPlot("Pseudocolor", "VITESSE_SOM_X_SOM_dom")
AddPlot("Pseudocolor", "VITESSE_X_FACES_dom_dual")

p.minFlag = 1
p.min = 0.9
p.maxFlag = 1
p.max = 1.1
p.centering = p.Natural
p.colorTableName = "hot_desaturated"
p.invertColorTable = 0
p.legendFlag = 0
p.lightingFlag = 0
SetPlotOptions(p)

AddOperator("Slice", 0)
s.originType = s.Percent
s.originPercent = 0
s.axisType = s.XAxis
s.project2d = 1
SetOperatorOptions(s, 0)

DrawPlots()

# Parametres de la legende
#-------------------------

nom_champ = GetPlotList().GetPlots(0).plotName
legende = GetAnnotationObject(nom_champ)
legende.active = 1
legende.managePosition = 0
legende.position = (0.1, 0.85)
legende.xScale = 1
legende.yScale = 0.6
legende.numberFormat = "%# -9.2f"
legende.fontHeight = 0.04
legende.drawTitle = 0
legende.drawMinMax = 0
legende.orientation = legende.VerticalRight  # VerticalRight, VerticalLeft, HorizontalTop, HorizontalBottom
legende.controlTicks = 1
legende.numTicks = 3
legende.fontFamily = legende.Times

# Parametres de la vue
#---------------------

v.windowCoords = (0, 0.1, 0, 0.1)
v.viewportCoords = (0.1, 0.9, 0.1, 0.9)
SetView2D(v)

# Parametres d'enregistrement d'images
#-------------------------------------

#InvertBackgroundColor()

sw = SaveWindowAttributes()
sw.outputToCurrentDirectory = 1
#sw.outputDirectory = "IMAGES"
sw.fileName = "image"
sw.format = sw.PNG
sw.width = 1024
sw.resConstraint = sw.ScreenProportions
SetSaveWindowAttributes(sw)

SaveWindow()

exit()
