execfile("instructions.py")

OpenDatabase("./test_convection_masse_lata_1.sauv.lata", 0)

AddPlot("Pseudocolor", "TEMPERATURE_ELEM_grid_geom2", 1, 0)
ChangeActivePlotsVar("__EXPRESSION_TO_PLOT__")
DrawPlots()
Query("Weighted Variable Sum")
moyenne_d_rho=GetQueryOutputValue()

f1=open('./integrale_d_rho.txt', 'w')
f1.write("integrale de d_rho_ :  " + str(moyenne_d_rho) + "\n")
f1.close

exit()
