execfile("instructions.py")

OpenDatabase("./test_diffusivite_turbulente_lata_1.sauv.lata", 0)

AddPlot("Curve", "operators/Lineout/TEMPERATURE_ELEM_grid_geom2", 1, 0)
ChangeActivePlotsVar("__EXPRESSION_TO_PLOT__")
LineoutAtts = LineoutAttributes()
LineoutAtts.point1 = (0.5, 0.5, 0.)
LineoutAtts.point2 = (0.5, 0.5, 1.)
LineoutAtts.interactive = 0
LineoutAtts.ignoreGlobal = 0
LineoutAtts.samplingOn = 0
LineoutAtts.numberOfSamplePoints = 50
LineoutAtts.reflineLabels = 0
SetOperatorOptions(LineoutAtts, 0)
DrawPlots()

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "./"
SaveWindowAtts.fileName = "__EXPRESSION_TO_PLOT__"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.CURVE  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 352
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

exit()
