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


lo=dir()
# Valeurs par defaut : 
nb_ts = 1       # frequence des pas de temps tires : Utilise dans range(1,TimeSliderGetNStates(), nb_ts)
ts_im = 1       # timestep chosen to buid images.
ts_start = 0    # timestep at which the loop starts
ts_end   = 10000# timestep at which the loop ends
p_filtre="p7"   # Choix de la taille du filtre
name_case="TCR" # Choix de la base de donnees
subname_case="Ja-1" # Choix du nombre de jacob du cas etudie.
# Entree des parametres physiques necessaires.
nu_l=2.e-5      # nu du liquide
k_l=0.2684      # conductivite du liquide
L_vap=-1.0e3    # Chaleur latente d'evaporation
T_sat= 0.       # Temperature de saturation
T_inf= 1.       # Temperature de bulk
Vinlet= 0.      # Valeur de la vitesse du referentiel si on travaille dans un referentiel mobile.
PAR=1           # Calcul paralelle : 0=non, 1=oui
GRADT=0         # EST-ce que le gradT est postTT ou seulement le mp : 0=non, 1=oui : Permet de savoir la variable e utiliser pour creer Grad_T_fois_surf_elem
DIR="./"
#DIR="/work/gb218285/Calculs/v154-pch/"+name_case+"/mmWe7.5-Reb150-"+subname_case+"-Pr1.474-n400-nitd25-dt1e-6-RelaxBary0.06-liss0-pas1.e-4/"
DIM="2D_axi"

print argv

# Quelques variables globales :

print "variables modifiables :" ,
for k in  dir():
   if k not in lo:
      print k,
      pass
   pass
print

for i in range(len(argv)):
   if argv[i]=="-val":
      exec(argv[i+1])
      i=i+1
      pass
   pass

print "ts_start = ",ts_start
print "ts_end = ",ts_end
print "nb_ts = ",nb_ts
print "ts_im = ",ts_im
print "p_filtre = ", p_filtre
print "name_case = ", name_case
print "subname_case = ", subname_case
print "DIR = ", DIR
print "nu_l = ", nu_l
print "k_l = ", k_l
print "L_vap = ", L_vap
print "T_sat = ", T_sat
print "T_inf = ", T_inf
print "PAR = ", PAR
print "GRADT = ", GRADT
print "Vinlet = ", Vinlet
print "DIM = ", DIM
if (DIM != "2D_axi") and (DIM != "3D"):
   raise Exception("variable DIM is not matching a predefined value")
   
#PathDats =DIR+"PostTT/"+p_filtre+"/Datas/"
PathDats =DIR

# On importe la bibliotheque mathematique pour utiliser le nombre Pi
import math

# On importe la base de donnees dans visit
print("Ready to import database...")
db=""
if os.path.isfile(DIR+"latas/compil.lata"):
   db=DIR+"latas/compil.lata"
elif os.path.isfile(DIR+"lata/post.lata"):
   db=DIR+"lata/post.lata"
else:
   import glob
   ldb=glob.glob(DIR+"*lata")
   if (len(ldb) == 0): raise Exception("Database missing")
   db=ldb[0]

OpenDatabase(db)
print("database %s imported..."%db)

plots = []
if (ts_end>TimeSliderGetNStates()):
   ts_end = TimeSliderGetNStates()
   pass

list_ts = range(ts_start,int(ts_end),nb_ts) # On centralise ici la liste des timestep e traiter.
print list_ts
pi = math.pi

DefineScalarExpression("xinterf", "coord(INTERFACES)[0]")
DefineScalarExpression("yinterf", "coord(INTERFACES)[1]")
if (DIM == "3D"):
   DefineScalarExpression("zinterf", "coord(INTERFACES)[2]")
   DefineScalarExpression("sinterf", "area(INTERFACES)")
   DefineScalarExpression("vol", "volume(dom)")
   DefineScalarExpression("vol", "volume(dom)")
else:
   DefineScalarExpression("sinterf", "min_edge_length(INTERFACES)*%g*2*coord(INTERFACES)[0]"%(pi))
   DefineScalarExpression("vol", "VOLUME_MAILLE_ELEM_dom")
   pass

# On definit une nouvelle variable (chi_v est l'indicatrice de la phase vapeur)
DefineScalarExpression("chi_v","1-INDICATRICE_INTERF_ELEM_dom")
DefineScalarExpression("vol_vap", "chi_v*vol")
DefineScalarExpression("volvitesse_x","VITESSE_X_ELEM_dom*vol_vap")
DefineScalarExpression("volvitesse_y","VITESSE_Y_ELEM_dom*vol_vap")
DefineScalarExpression("volvitesse_z","VITESSE_Z_ELEM_dom*vol_vap")
# On construit le (gradT.n) soit e partir de mpoint, soit e partir de gradT selon ce dont on dispose.
if GRADT==1 :
   DefineScalarExpression("Grad_T_fois_surf_elem", "pos_cmfe(<[0]id:GRAD_TEMP_THERMIQUE_ELEM_dom>,INTERFACES,-100.)*sinterf")
else :
   # DefineScalarExpression("Grad_T_fois_surf_elem",str(L_vap)+"/"+str(k_l)+"*pos_cmfe(<[0]id:MPOINT_THERMIQUE_ELEM_dom>,INTERFACES,-100.)*sinterf")
   DefineScalarExpression("Grad_T_fois_surf_elem",str(L_vap)+"/"+str(k_l)+"*pos_cmfe(<[0]id:MPOINT_THERMIQUE_ELEM_dom>,INTERFACES,-100.)*sinterf")

if PAR==1 :
   DefineScalarExpression("facettes_locales", "eq(PE_LOCAL_ELEM_INTERFACES, PE_ELEM_INTERFACES)") # Pour eliminer facettes virtuelles

# L'idee c'est d'avoir 2 plots : 
#    0 -> pour le maillage eulerien
#    1 -> pour le maillage lagrangien

# On cree un trace type pseudocolor
AddPlot("Pseudocolor","COMPO_CONNEXE_ELEM_INTERFACES")

# On affiche l'interface
#AddPlot("Mesh","INTERFACES")
#AddPlot("Pseudocolor","INTERFACES") # Mais sous la forme d'un pseudocolor (donc les triangles sont pleins et colores... )
AddPlot("Pseudocolor","facettes_locales")
res=[]
if PAR==1 :
   AddOperator("Threshold")
   ThresholdAtts = ThresholdAttributes()
   #ThresholdAtts.outputMeshType = 0
   ThresholdAtts.listedVarNames = ("facettes_locales", "COMPO_CONNEXE_ELEM_INTERFACES", "yinterf")
   ThresholdAtts.zonePortions = (1, 1, 1)
   ThresholdAtts.lowerBounds = (1, -1e+37, -1e+37)
   ThresholdAtts.upperBounds = (1e+37, 1e+37, 1e+37)
   ThresholdAtts.defaultVarName = "facettes_locales"
   ThresholdAtts.defaultVarIsScalar = 1
   SetOperatorOptions(ThresholdAtts)
   pass

AddOperator("Box", 1)
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.Some  # Some, All
BoxAtts.minx = -1
BoxAtts.maxx = 1
BoxAtts.miny = -1
BoxAtts.maxy = 1
BoxAtts.minz = 0
BoxAtts.maxz = 2
BoxAtts.inverse = 0
#  This function takes an operator attributes object as the first argument. The second argument, which is optional is the active operator index. When it is specified, the operator attributes will only be applied to specified active operator. The third argument, which is also optional, is a flag that tells VisIt if the operator attributes should be applied to all plots instead of just the selected plots.
SetOperatorOptions(BoxAtts, 1,1) # To all

DrawPlots()
ListPlots()
plots = plots + [0] + [1]

global D_MAX
D_MAX = 1e+37
def mod_ThresholdAtts(ThresholdAtts, compo_min=-D_MAX, compo_max=D_MAX, y_min=-D_MAX, y_max=D_MAX):
   #ThresholdAtts.outputMeshType = 0
   ThresholdAtts.listedVarNames = ("default", "facettes_locales", "COMPO_CONNEXE_ELEM_INTERFACES", "yinterf")
   ThresholdAtts.zonePortions = (1, 1, 1, 1)
   ThresholdAtts.lowerBounds = (-1e+37, 1, compo_min, y_min)
   ThresholdAtts.upperBounds = (1e+37, 1e+37, compo_max, y_max)
   ThresholdAtts.defaultVarName = "INTERFACES"
   ThresholdAtts.defaultVarIsScalar = 0
   SetOperatorOptions(ThresholdAtts) # Probably just once
   pass

def mod_BoxAtts(BoxAtts, y_min=-D_MAX, y_max=D_MAX, z_min=-D_MAX, z_max=D_MAX):
   BoxAtts.amount = BoxAtts.Some  # Some, All
   BoxAtts.minx = -1
   BoxAtts.maxx = 1
   BoxAtts.miny = y_min
   BoxAtts.maxy = y_max
   BoxAtts.minz = z_min
   BoxAtts.maxz = z_max
   BoxAtts.inverse = 0
   SetOperatorOptions(BoxAtts, 1,1) # To all
   pass

# On cree le fichier contenant les resultats
nb_max_authorized_bubbles = 5
dy_tol = 10e-6 # Tolerance sur le domaine a garder... Au moins une taille de mailles en z pour etre sur de garder l'integralite des mailles euleriennes diphasiques (sinon, elles seront coupees) 
filename = "resultats"
ff = []
for i in range(nb_max_authorized_bubbles): 
   f = open(PathDats+filename+"_%d.dat"%(i), "w")
   ff.append(f)
   del f
   ff[i].write("# 1-index |  2-Time  |  3-S(t)  |  4-S(t)/S(0)  | 5-S(t)/S_sphere(t) |  6-V(t)  |  7-V(t)/V(t=0)  | 8-D_eq  | 9-Vit_vap_x  | 10-Vit_vap_y  | 11-Vit_vap_z  | 12-Norme_Vit_vap  | 13-Re_b=Norme_Vit_vap*D_eq/nu_l  | 14-Nusselt | 15-x_g | 16-y_g | 17-z_g  \n")
   pass

def reset():
      # Reset operators to initial state (utilise les valeurs par defaut des fonctions) : 
      SetActivePlots(1) # Activer le plot qui contient le Seuil. 
      mod_ThresholdAtts(ThresholdAttributes()) 
      mod_BoxAtts(BoxAttributes())
      pass

# On fait une boucle sur l'ensemble des pas de temps
# fonction range(a,b,c) : a = premiere valeur, b = derniere valeur, c = pas entre deux valeurs consecutives 
# for state in range(0,TimeSliderGetNStates(),1):
print('BEGIN TIMELOOP')
for state in list_ts: 
   #------------------------------------
   #         Pour la dynamique
   #------------------------------------
   reset() # On fait le reset avant de changer le timestep pour eviter de faire une "box" sur une region 
   #         ou il n'y a plus d'interface ce qui conduirait a la levee d'un warning. 
   SetTimeSliderState(state)
   Query("Time")
   time=GetQueryOutputValue()
   SetActivePlots(1)
   ChangeActivePlotsVar("COMPO_CONNEXE_ELEM_INTERFACES")
   Query("MinMax", use_actual_data=1)
   ib, jb =GetQueryOutputValue()
   ib,jb=int(ib),int(jb)
   nbulles = jb+1-ib
   lbulles =  range(ib,jb+1)
   print "*"*60, "\n", "At timestep ", state, " there is/are ", nbulles, " bubbles : ", lbulles
   if (nbulles > nb_max_authorized_bubbles): raise Exception("There are more bubbles than expected at most")
   for bulle in lbulles:
      reset()
      print "Dealing with bubble ", bulle,
      SetActivePlots(1)
      ChangeActivePlotsVar("yinterf")
      mod_ThresholdAtts(ThresholdAttributes(),compo_min=bulle, compo_max=bulle) 
      Query("MinMax", use_actual_data=1)
      ymin, ymax = GetQueryOutputValue()
      print " located between", ymin, " and ", ymax
      mod_BoxAtts(BoxAttributes(), y_min=ymin-dy_tol, y_max=ymax+dy_tol)
      SetActivePlots(0)
      ChangeActivePlotsVar("vol_vap")
      Query("Variable Sum")
      
      if state==list_ts[0] :
         volume_bulle_init=GetQueryOutputValue()
         volume_bulle=volume_bulle_init
      else: 
         volume_bulle=GetQueryOutputValue()
         pass
      # On calcul le diametre d'une sphere equivalente ayant le volume precedemment obtenu :
      # En 3D : V = 1/6*Pi*D^3
      diam_bulle_equiv = (volume_bulle*6./math.pi)**(1./3.)
      ChangeActivePlotsVar("volvitesse_x")
      Query("Variable Sum")
      vitesse_bulle_x=GetQueryOutputValue()/volume_bulle 
      ChangeActivePlotsVar("volvitesse_y")
      Query("Variable Sum")
      vitesse_bulle_y=GetQueryOutputValue()/volume_bulle
      vitesse_bulle_z=0
      if (DIM == "3D"):
         ChangeActivePlotsVar("volvitesse_z")
         Query("Variable Sum")
         # Suivant z, si on impose au liquide une vitesse Vinlet suivant z, alors la vrai difference de vitesse est :
         # dv = - V_referentiel + V_apparente
         vitesse_bulle_z=GetQueryOutputValue()/volume_bulle - Vinlet
      else:
         pass
      print "La vitesse de la bulle est de %g, %g ,%g" % (vitesse_bulle_x,vitesse_bulle_y,vitesse_bulle_z)
      # On en deduit la norme de la vitesse
      vitesse_bulle_norme=(vitesse_bulle_x**2.+vitesse_bulle_y**2+vitesse_bulle_z**2)**(1./2.)
      # On en deduit le nombre de Reynolds
      Re_bulle= vitesse_bulle_norme*diam_bulle_equiv/nu_l
      #------------------------------------
      #     Pour le centre de gravite
      #------------------------------------
      SetActivePlots(1)
      m=[]
      pos_bulle=[]
      ChangeActivePlotsVar("INTERFACES")
      #if PAR==1 :
      #   AddOperator("Threshold") 
      #   SetOperatorOptions(ThresholdAtts)
      Query("Centroid",1,0,"Default")
      m=GetQueryOutputValue()
      for toto in m:
         pos_bulle.append(toto)
      #if PAR==1 :
      #   RemoveLastOperator()
      print "Barycentre : ", m
      try:
         ChangeActivePlotsVar("sinterf")
         Query("Variable Sum")
      except:
         Query("3D surface area")
      if state==list_ts[0] :
         surface_init=GetQueryOutputValue()
         surface=surface_init
      else :
         surface=GetQueryOutputValue()
      print "La surface de la bulle est de %g " % surface
      # On calcul la surface d'une sphere equivalente ayant le volume precedement obtenu : 
      surface_equiv = (volume_bulle*6.)**(2./3.)*(math.pi)**(1./3.)
      print "Le volume (resp. surf) de la bulle est de %g m^3 (resp m^2), la surf equiv (resp. perimetre) est %g et le diametre equivalent est de %g m" % (volume_bulle,surface_equiv,diam_bulle_equiv)
      #------------------------------------
      #         Pour la thermique
      #------------------------------------
      ChangeActivePlotsVar("Grad_T_fois_surf_elem")
      #if PAR==1 :
      #   AddOperator("Threshold")
      Query("Variable Sum")
      Somme_Grad_T=GetQueryOutputValue()
      #if PAR==1 :
      #   RemoveLastOperator()
      # On en deduit le nombre de Nusselt
      # Dans le cas general, Nu = (gradT.n)_moyen * D_eq / (T_inf - T_sat)
      # On peut donc le calculer par : Nu = (gradT.n)*Aire *D_eq / (Aire * (T_inf - T_sat))
      # En 3D, Boris utilisait : Nu = (gradT.n)*Aire/(Pi*D*(T_inf - T_sat)), car pour une sphere Aire = Pi*D^2.
      Nusselt=Somme_Grad_T*diam_bulle_equiv/(surface_equiv*(T_inf-T_sat))
      #------------------------------------
      # Ecriture des resultats dans un fichier
      #------------------------------------
      s="%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n" %(state, \
                                         time, \
                                         surface, \
                                         surface/surface_init, \
                                         surface/surface_equiv, \
                                         volume_bulle, \
                                         volume_bulle/volume_bulle_init, \
                                         diam_bulle_equiv, \
                                         vitesse_bulle_x, \
                                         vitesse_bulle_y, \
                                         vitesse_bulle_z, \
                                         vitesse_bulle_norme, \
                                         Re_bulle, \
                                         Nusselt, \
                                         pos_bulle[0], \
                                         pos_bulle[1], \
                                         pos_bulle[2], \
                                         )
      ff[bulle].write(s)
      pass
   pass
#f.close()

for i in range(nb_max_authorized_bubbles): 
   ff[i].close()
   pass

print('THE END')
exit(0)
# T_inf= float(raw_input("Valeur de Tinfini ?"))

