#!$TRUST_ROOT/exec/VisIt/bin/visit -nowin -cli -s 
OpenDatabase("./lata/post.lata", 0)
DefineScalarExpression("dTdni", "recenter(pos_cmfe(<[0]id:TEMPERATURE_GRAD_THERMIQUE_ELEM_dom>, INTERFACES, -1e8), \"nodal\")")
AddPlot("Pseudocolor", "dTdni", 1, 1)
DrawPlots()
#DefineScalarExpression("xi", "coord(INTERFACES)[0]")
#DefineScalarExpression("yi", "coord(INTERFACES)[1]")
nb = TimeSliderGetNStates()
SetTimeSliderState(nb-1)
ExportDBAtts = ExportDBAttributes()
ExportDBAtts.allTimes = 0
ExportDBAtts.db_type = "XYZ"
ExportDBAtts.db_type_fullname = "XYZ_1.0"
ExportDBAtts.filename = "dTdni"
ExportDBAtts.dirname = "."
ExportDBAtts.variables = ("dTdni")
ExportDBAtts.writeUsingGroups = 0
ExportDBAtts.groupSize = 48
ExportDBAtts.opts.types = ()
ExportDBAtts.opts.help = ""
ExportDatabase(ExportDBAtts)
exit()
