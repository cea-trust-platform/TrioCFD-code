#visit -cli -s Visus.py -nowin

Maillages =     [
                'mesh_hexa_1',
                'mesh_hexa_2',
                'mesh_hexa_3',
                'mesh_hexa_4',
                'mesh_hexa_5',
                'mesh_tetra_0',
                'mesh_tetra_1',
                'mesh_tetra_2',
                'mesh_tetra_3',
                'mesh_tetra_4',
                'mesh_tetra_5',
                'mesh_tetra_6'
                ]


N=len(Maillages)

a = AnnotationAttributes()
v = GetView3D()
sw = SaveWindowAttributes()

# Parametres du rendu
#--------------------

r = RenderingAttributes()
r.antialiasing = 1
SetRenderingAttributes(r)

# Parametres des annotations
#---------------------------

a.axes2D.visible = 1
a.axes2D.autoSetTicks = 1
a.axes2D.autoSetScaling = 1
a.axes2D.lineWidth = 0
a.axes2D.xAxis.title.visible = 0
a.axes2D.xAxis.title.userUnits = 1
a.axes2D.xAxis.label.visible = 1
a.axes2D.xAxis.tickMarks.visible = 1
a.axes2D.xAxis.tickMarks.majorMinimum = 0
a.axes2D.xAxis.tickMarks.majorMaximum = 1
a.axes2D.xAxis.tickMarks.minorSpacing = 1
a.axes2D.xAxis.tickMarks.majorSpacing = 1
a.axes2D.yAxis.title.visible = 0
a.axes2D.yAxis.title.userUnits = 1
a.axes2D.yAxis.label.visible = 1
a.axes2D.yAxis.tickMarks.visible = 1
a.axes2D.yAxis.tickMarks.majorMinimum = 0
a.axes2D.yAxis.tickMarks.majorMaximum = 1
a.axes2D.yAxis.tickMarks.minorSpacing = 1
a.axes2D.yAxis.tickMarks.majorSpacing = 1
a.userInfoFlag = 0
a.databaseInfoFlag = 1
a.timeInfoFlag = 0
a.databaseInfoTimeScale = 1
a.databaseInfoTimeOffset = 0
a.legendInfoFlag = 0
SetAnnotationAttributes(a)

# Parametres de la vue
#---------------------

v.viewNormal = (-0.56819, 0.446318, 0.691347)
v.focus = (0.5, 0.5, 0.5)
v.viewUp = (0.252497, 0.894169, -0.369739)
v.viewAngle = 30
v.parallelScale = 0.866025
v.nearPlane = -1.73205
v.farPlane = 1.73205
v.imagePan = (0, 0)
v.imageZoom = 0.885164
v.perspective = 1
v.eyeAngle = 2
v.centerOfRotationSet = 0
v.centerOfRotation = (0.5, 0.5, 0.5)
v.axis3DScaleFlag = 0
v.axis3DScales = (1, 1, 1)
v.shear = (0, 0, 1)
v.windowValid = 1
SetView3D(v)

# Parametres d'enregistrement d'images
#-------------------------------------

sw.outputToCurrentDirectory = 1
sw.family = 1
sw.format = sw.PNG
sw.width = 1024
sw.height = 1024
sw.screenCapture = 0
sw.saveTiled = 0
sw.quality = 80
sw.progressive = 0
sw.binary = 0
sw.stereo = 0
sw.compression = sw.PackBits
sw.forceMerge = 0
sw.resConstraint = sw.ScreenProportions
sw.advancedMultiWindowSave = 0
SetSaveWindowAttributes(sw)

# Ajout des maillages
#--------------------

for j in range(N):
    dom=Maillages[j]
    maillage=dom+'.med'
    OpenDatabase(maillage)
    AddPlot("Mesh",dom)
    DrawPlots()
    sw.fileName = dom
    SetSaveWindowAttributes(sw)
    SaveWindow()
    DeleteActivePlots()

exit()
