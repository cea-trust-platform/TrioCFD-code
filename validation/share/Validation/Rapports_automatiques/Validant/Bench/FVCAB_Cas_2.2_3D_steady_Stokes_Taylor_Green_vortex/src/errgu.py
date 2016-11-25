#visit -cli -s errgu.py -nowin

Maillages =	[
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

for i in range(N):
	maillage=Maillages[i]
	db=maillage+'/Cas.lata'
	fic=maillage+'/val.dat'
	OpenDatabase(db)
	SetTimeSliderState(TimeSliderGetNStates() - 1)
	DefineVectorExpression("uex", "{-2*cos(2*3.1415926536*coord(dom)[0])*sin(2*3.1415926536*coord(dom)[1])*sin(2*3.1415926536*coord(dom)[2]), sin(2*3.1415926536*coord(dom)[0])*cos(2*3.1415926536*coord(dom)[1])*sin(2*3.1415926536*coord(dom)[2]), sin(2*3.1415926536*coord(dom)[0])*sin(2*3.1415926536*coord(dom)[1])*cos(2*3.1415926536*coord(dom)[2])}")
	DefineScalarExpression("uex_1", "uex[0]")
	DefineScalarExpression("uex_2", "uex[1]")
	DefineScalarExpression("uex_3", "uex[2]")
	DefineScalarExpression("norme_gradient_uex_carre", "gradient(uex_1)[0]^2+gradient(uex_1)[1]^2+gradient(uex_1)[2]^2+gradient(uex_2)[0]^2+gradient(uex_2)[1]^2+gradient(uex_2)[2]^2+gradient(uex_3)[0]^2+gradient(uex_3)[1]^2+gradient(uex_3)[2]^2")
	DefineVectorExpression("u_uex", "{VITESSE_X_SOM_dom-uex_1,VITESSE_Y_SOM_dom-uex_2,VITESSE_Z_SOM_dom-uex_3}")
	DefineScalarExpression("u_uex_1", "u_uex[0]")
	DefineScalarExpression("u_uex_2", "u_uex[1]")
	DefineScalarExpression("u_uex_3", "u_uex[2]")
	DefineScalarExpression("norme_gradient_u_uex_carre", "gradient(u_uex_1)[0]^2+gradient(u_uex_1)[1]^2+gradient(u_uex_1)[2]^2+gradient(u_uex_2)[0]^2+gradient(u_uex_2)[1]^2+gradient(u_uex_2)[2]^2+gradient(u_uex_3)[0]^2+gradient(u_uex_3)[1]^2+gradient(u_uex_3)[2]^2")
	AddPlot("Pseudocolor","norme_gradient_uex_carre")
	DrawPlots()
	Query("Average Value")
	val1 = GetQueryOutputValue()
	DeleteActivePlots()
	AddPlot("Pseudocolor","norme_gradient_u_uex_carre")
	DrawPlots()
	Query("Average Value")
	val2 = GetQueryOutputValue()
	DeleteActivePlots()
	filout = open(fic,'w')
	filout.write(str(val1)+'\t')
	filout.write(str(val2))
	filout.close()

exit()



