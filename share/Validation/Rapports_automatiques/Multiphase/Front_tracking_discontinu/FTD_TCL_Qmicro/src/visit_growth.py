# -*- coding:Utf-8 -*-
# Ligne precedente permet l'utilisation des caracteres accentues
"""
      Script Python permettant de calculer :
                            - vitesse d'une bulle de vapeur
                            - reynolds de bulle
                            - nusselt de bulle
                            - centre de gravite des bulles (avec la fonction "CENTROID")
 Teste avec VisIt 1.7

"""
# ----------------------------------------------------------------------
# Script Python permettant de calculer des correlations (Teste avec VisIt 1.8.1)
# -----------------------------------------------------------------------
from sys import *

ts_start = 1 # On zappe le 0
nb_ts = 1
ts_end = 100000
PAR=0
# On importe la bibliotheque mathematique pour utiliser le nombre Pi
import math
pi = math.pi

# On importe la base de donnees dans visit
print("Ready to import database...")
db=""
if os.path.isfile("lata/post.lata"):
   db="lata/post.lata"
else:
   import glob
   ldb=glob.glob("*/*lata")
   if (len(ldb) == 0): raise Exception("Database missing")
   db=ldb[0]

OpenDatabase(db)
print("database %s imported..."%db)

plots = []
if (ts_end>TimeSliderGetNStates()):
   ts_end = TimeSliderGetNStates()
   pass

list_ts = range(ts_start,int(ts_end),nb_ts) # On centralise ici la liste des timestep e traiter.
print(list_ts)

DefineScalarExpression("xinterf", "coord(INTERFACES)[0]")
DefineScalarExpression("yinterf", "coord(INTERFACES)[1]")

DefineScalarExpression("sinterf", "min_edge_length(INTERFACES)*%g*2*coord(INTERFACES)[0]"%(pi))
DefineScalarExpression("vol", "VOLUME_MAILLE_ELEM_dom")

# On definit une nouvelle variable (chi_v est l'indicatrice de la phase vapeur)
DefineScalarExpression("chi_v","1-INDICATRICE_INTERF_ELEM_dom")
DefineScalarExpression("vol_vap", "chi_v*vol")
DefineScalarExpression("volvitesse_x","VITESSE_X_ELEM_dom*vol_vap")
DefineScalarExpression("volvitesse_y","VITESSE_Y_ELEM_dom*vol_vap")
DefineScalarExpression("volvitesse_z","VITESSE_Z_ELEM_dom*vol_vap")
# On construit le (gradT.n) soit e partir de mpoint, soit e partir de gradT selon ce dont on dispose.
DefineScalarExpression("mpai_ELEM_dom", "MPOINT_THERMIQUE_ELEM_dom*INTERFACIAL_AREA_ELEM_dom")

if PAR==1 :
   DefineScalarExpression("facettes_locales", "eq(PE_LOCAL_ELEM_INTERFACES, PE_ELEM_INTERFACES)") # Pour eliminer facettes virtuelles

# L'idee c'est d'avoir 2 plots : 
#    0 -> pour le maillage eulerien
#    1 -> pour le maillage lagrangien

# On cree un trace type pseudocolor
AddPlot("Pseudocolor","mpai_ELEM_dom")

# On affiche l'interface
#AddPlot("Mesh","INTERFACES")
#AddPlot("Pseudocolor","INTERFACES") # Mais sous la forme d'un pseudocolor (donc les triangles sont pleins et colores... )
AddPlot("Pseudocolor","xinterf")
res=[]
if 0 :
   AddOperator("Threshold")
   ThresholdAtts = ThresholdAttributes()
   #ThresholdAtts.outputMeshType = 0
   ThresholdAtts.listedVarNames = ("facettes_locales", "COMPO_CONNEXE_ELEM_INTERFACES", "yinterf")
   ThresholdAtts.zonePortions = (1, 1, 1)
   ThresholdAtts.lowerBounds = (1, -1e+37, -1e+37)
   ThresholdAtts.upperBounds = (1e+37, 1e+37, 20e-6)
   ThresholdAtts.defaultVarName = "facettes_locales"
   ThresholdAtts.defaultVarIsScalar = 1
   SetOperatorOptions(ThresholdAtts)
   pass

DrawPlots()
SetActivePlots(1)
AddOperator("Box", 0)
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.Some  # Some, All
BoxAtts.minx = 0
BoxAtts.maxx = 1
BoxAtts.miny = 0
BoxAtts.maxy = 2e-05
BoxAtts.minz = 0
BoxAtts.maxz = 1
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 0, 0)


DrawPlots()
ListPlots()
plots = plots + [0] + [1]

global D_MAX
D_MAX = 1e+37

# On fait une boucle sur l'ensemble des pas de temps
# fonction range(a,b,c) : a = premiere valeur, b = derniere valeur, c = pas entre deux valeurs consecutives 
# for state in range(0,TimeSliderGetNStates(),1):
print('BEGIN TIMELOOP')
f = open("results.dat", 'w')
f.write("# State time xcl mpai Vb dVbdt_from_mpai dVbdt_ana\n")

Qcl = 50.
Lvap = -2.256e6 # J/kg
rhov=800.  # kg/m3
for state in list_ts: 
   #------------------------------------
   SetTimeSliderState(state)
   Query("Time")
   time=GetQueryOutputValue()
   SetActivePlots(0)
   ChangeActivePlotsVar("mpai_ELEM_dom")
   Query("Variable Sum")
   mpai=GetQueryOutputValue()
   
   ChangeActivePlotsVar("vol_vap")
   Query("Variable Sum")
   Vb=GetQueryOutputValue()
   
   SetActivePlots(1)
   ChangeActivePlotsVar("xinterf")
   Query("Min", use_actual_data=1) # 1 -> to get the min only after the box
   xcl=GetQueryOutputValue()
   
   dVbdt_ana = Qcl* 2.*math.pi *xcl / (rhov*Lvap)# m3/(s)
   dVbdt_from_mpai = mpai/rhov # m3/(s)
   print("timestep ", state, " time ", time, " xcl ", xcl, " integral(mpai) ", mpai, " Vb ", Vb, " dVbdt ", dVbdt_from_mpai, " dVbdt_ana ", dVbdt_ana)
   #------------------------------------
   # Ecriture des resultats dans un fichier
   #------------------------------------
   s="%d %g %g %g %g %g %g\n" %(state, \
                                         time, \
                                         xcl, \
                                         mpai, \
                                         Vb, \
                                         dVbdt_from_mpai, \
                                         dVbdt_ana \
                                         )
   f.write(s)
   pass
f.close()

print('THE END')
exit(0)
# T_inf= float(raw_input("Valeur de Tinfini ?"))
