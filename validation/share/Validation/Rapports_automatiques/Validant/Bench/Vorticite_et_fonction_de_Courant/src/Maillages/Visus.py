#visit -cli -s Visus.py -nowin

Maillages =	[
		'mesh_cart_1',
		'mesh_cart_2',
		'mesh_cart_3',
		'mesh_cart_4',
		'mesh_cart_5',
		'mesh_cart_6',
		'mesh_cart_7',
		'mesh_quad_1',
		'mesh_quad_2',
		'mesh_quad_3',
		'mesh_quad_4',
		'mesh_quad_5',
		'mesh_quad_6',
		'mesh_quad_7',
		'mesh_tri_1',
		'mesh_tri_2',
		'mesh_tri_3',
		'mesh_tri_4',
		'mesh_tri_5',
		'mesh_tri_6'
		]


N=len(Maillages)

a = AnnotationAttributes()
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
