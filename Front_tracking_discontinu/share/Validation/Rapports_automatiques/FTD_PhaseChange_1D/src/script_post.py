# -*- coding:Utf-8 -*-
"""
      Script Python permettant de calculer :
"""
#import numpy as np
import sys, re
from sys import *
import math
import visit

def getValueFromFic_triou(fic, key):
   f = open(fic, 'r')
   lines = f.readlines()
   f.close()
   nb=len(lines)
   rc=re.compile(key+"[\s]*(?P<value>[\-]?[\d]*.[\d]*[eE]?[+\-]?[\d]*)")
   #rc=re.compile(key+"[\s]*(?P<value>[\-]?[\d]*.[\d]*[eE]?[+\-]?[\d]*)")
   #print fic, key, delta
   for i, st in enumerate(lines):
      m=rc.match(st.strip())
      if m:
         val = float(m.group("value"))
         break
      pass
   if i == nb- 1 :
      print("On a rien trouvee??")
      raise Exception("Etonnant, non?")
      pass
   return val


if (1):	
	print("gogogo...", os.getcwd())
	Sur=144.e-10
	dt_sauv=9.e-8
	dt2=1.e-8
	aa = getValueFromFic_triou("adiab/adia_remesh.data", "dt_max")
	print(aa==dt2)
	f=100000000.
	rhov=1000.
	cpv=100.
	kv=0.4
	gradTv=0.
	rhol=1000.
	cpl=600.
	Lv_l=3.
	kl=0.6
	gradTl=-2.5e5
	flux_ent=15.e4
	L_maille=1.e-5

	dt=dt_sauv
	print(("Sur=",Sur))
	print(("dt=",dt))
	print(("dt2=",dt2))
	print(("f=",f))
	print(("rhov=",rhov))
	print(("cpv=",cpv))
	print(("kv=",kv))
	print(("gradTv=",gradTv))
	print(("rhol=",rhol))
	print(("cpl=",cpl))
	print(("Lv_l=",Lv_l))
	print(("kl=",kl))
	print(("gradTl=",gradTl))
	print(("flux_ent=",flux_ent))
	print(("L_maille=",L_maille))

	OpenDatabase("adiab/lata/post.lata")

	ts_start=0
	ts_end = TimeSliderGetNStates()
	nb_ts=1
	list_ts = list(range(ts_start,ts_end,nb_ts)) 

	DefineScalarExpression("Tv_ELEM_dom","TEMPERATURE_THERMIQUE_VAPEUR_ELEM_dom*(1-INDICATRICE_INTERF_ELEM_dom)")
	#DefineScalarExpression("Tl_ELEM_dom","TEMPERATURE_THERMIQUE_ELEM_dom*INDICATRICE_INTERF_ELEM_dom")
	DefineScalarExpression("MPOINT_ELEM_INTERFACES","pos_cmfe(<[0]id:TEMPERATURE_MPOINT_ELEM_dom>,INTERFACES,0.)")
	DefineScalarExpression("MPOINTV_ELEM_INTERFACES","pos_cmfe(<[0]id:TEMPERATURE_MPOINTV_ELEM_dom>,INTERFACES,0.)")
	'''DefineScalarExpression("Pressure_Power","GRADIENT_PRESSION_X_ELEM_dom*VITESSE_X_ELEM_dom")'''
	#DefineScalarExpression("CHALEUR_MASSIQUE_ELEM_dom","cp*rho*INDICATRICE_INTERF_ELEM_dom*TEMPERATURE_THERMIQUE_ELEM_dom")
	DefineScalarExpression("Tl_ELEM_dom","if(or(INDICATRICE_INTERF_ELEM_dom=1, INDICATRICE_INTERF_ELEM_dom=0), TEMPERATURE_THERMIQUE_ELEM_dom*INDICATRICE_INTERF_ELEM_dom, -(INDICATRICE_INTERF_ELEM_dom)^2/2*TEMPERATURE_GRAD_THERMIQUE_ELEM_dom*%g)"%(L_maille))	

	# Temps initial : 
	SetTimeSliderState(0)
	AddPlot("Mesh","INTERFACES")
	DrawPlots()
	Query("Centroid")
	x0,y0,z0=GetQueryOutputValue()
	DeleteActivePlots()
	AddPlot("Pseudocolor","Tv_ELEM_dom",1,1)
	DrawPlots()
	Query("Weighted Variable Sum")
	Tv0=GetQueryOutputValue()
	ChangeActivePlotsVar("Tl_ELEM_dom")
	Query("Weighted Variable Sum")
	Tl0=GetQueryOutputValue()


	# Temps final 
	SetTimeSliderState(9)
	DeleteActivePlots()
	AddPlot("Mesh","INTERFACES")
	DrawPlots()
	Query("Centroid")
	x1,y1,z1=GetQueryOutputValue()
	DeleteActivePlots()
	AddPlot("Pseudocolor","MPOINT_ELEM_INTERFACES",1,1)
	DrawPlots()
	Query("Weighted Variable Sum")
	mpl1=GetQueryOutputValue()
	ChangeActivePlotsVar("MPOINTV_ELEM_INTERFACES")
	Query("Weighted Variable Sum")
	mpv1=GetQueryOutputValue()
	ChangeActivePlotsVar("Tv_ELEM_dom")
	Query("Weighted Variable Sum")
	Tv1=GetQueryOutputValue()
	ChangeActivePlotsVar("Tl_ELEM_dom")
	Query("Weighted Variable Sum")
	Tl1=GetQueryOutputValue()
	'''ChangeActivePlotsVar("Pressure_Power")
	Query("Weighted Variable Sum")
	PP=GetQueryOutputValue()'''
	ChangeActivePlotsVar("VITESSE_X_ELEM_dom")
	Query("Average value")
	U_cine=GetQueryOutputValue()
	print(((U_cine**2)*200.e-6*Sur/2*(rhol+rhov)/2))

	ChangeActivePlotsVar("Tv_ELEM_dom")
	AddOperator("Slice", 1)
	SliceAtts = SliceAttributes()
	SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
	SliceAtts.originPoint = (1.9999e-4, 0, 0)
	SliceAtts.originIntercept = 0
	SliceAtts.originPercent = 0
	SliceAtts.originZone = 0
	SliceAtts.originNode = 0
	SliceAtts.normal = (-1, 0, 0)
	SliceAtts.axisType = SliceAtts.XAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.upAxis = (0, 1, 0)
	SliceAtts.project2d = 1
	SliceAtts.interactive = 1
	SliceAtts.flip = 0
	SliceAtts.originZoneDomain = 0
	SliceAtts.originNodeDomain = 0
	SliceAtts.meshName = "INTERFACES"
	SliceAtts.theta = 90
	SliceAtts.phi = 0
	SetOperatorOptions(SliceAtts, 1)
	DrawPlots()
	Query("Average Value")
	Tvs1=GetQueryOutputValue()
	ChangeActivePlotsVar("VITESSE_X_ELEM_dom")
	Query("Average value")
	U_sor=-GetQueryOutputValue()
	'''ChangeActivePlotsVar("PRESSION_ELEM_dom")
	Query("Average value")
	P_sor=-GetQueryOutputValue()'''
	RemoveAllOperators(1)


	ChangeActivePlotsVar("Tl_ELEM_dom")
	AddOperator("Slice", 1)
	SliceAtts = SliceAttributes()
	SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
	SliceAtts.originPoint = (0.0001e-4, 0, 0)
	SliceAtts.originIntercept = 0
	SliceAtts.originPercent = 0
	SliceAtts.originZone = 0
	SliceAtts.originNode = 0
	SliceAtts.normal = (-1, 0, 0)
	SliceAtts.axisType = SliceAtts.XAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.upAxis = (0, 1, 0)
	SliceAtts.project2d = 1
	SliceAtts.interactive = 1
	SliceAtts.flip = 0
	SliceAtts.originZoneDomain = 0
	SliceAtts.originNodeDomain = 0
	SliceAtts.meshName = "INTERFACES"
	SliceAtts.theta = 90
	SliceAtts.phi = 0
	SetOperatorOptions(SliceAtts, 1)
	DrawPlots()
	Query("Average Value")
	Tle1=GetQueryOutputValue()
	ChangeActivePlotsVar("VITESSE_X_ELEM_dom")
	Query("Average value")
	U_ent=-GetQueryOutputValue()
	'''ChangeActivePlotsVar("PRESSION_ELEM_dom")
	Query("Average value")
	P_sor=-GetQueryOutputValue()'''
	RemoveAllOperators(1)

	f = open("adiab/adia_remesh_pb_Diffusion_chaleur.out")
	lines=f.readlines()
	f.close()
	st = lines[-1]
	l=st.split()
	tps=float(l[0])
	diff_sortie=float(l[1])
	Pui1=diff_sortie

	f = open("adiab/adia_remesh_pb_Convection_chaleur.out")
	lines=f.readlines()
	f.close()
	st = lines[-1]
	l=st.split()
	tps=float(l[0])
	conv_sortie=float(l[3])

	inter=(kv*gradTv-kl*gradTl)/Lv_l

	#Calcul de dx c