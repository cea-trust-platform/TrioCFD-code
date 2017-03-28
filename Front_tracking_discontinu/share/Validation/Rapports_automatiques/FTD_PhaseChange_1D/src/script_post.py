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
      print "On a rien trouvee??"
      raise Exception("Etonnant, non?")
      pass
   return val


if (1):	
	print "gogogo...", os.getcwd()
	Sur=144.e-10
	dt_sauv=9.e-8
	dt2=1.e-8
	aa = getValueFromFic_triou("adiab/adia_remesh.data", "dt_max")
	print aa==dt2
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
	print("Sur=",Sur)
	print("dt=",dt)
	print("dt2=",dt2)
	print("f=",f)
	print("rhov=",rhov)
	print("cpv=",cpv)
	print("kv=",kv)
	print("gradTv=",gradTv)
	print("rhol=",rhol)
	print("cpl=",cpl)
	print("Lv_l=",Lv_l)
	print("kl=",kl)
	print("gradTl=",gradTl)
	print("flux_ent=",flux_ent)
	print("L_maille=",L_maille)

	OpenDatabase("adiab/lata/post.lata")

	ts_start=0
	ts_end = TimeSliderGetNStates()
	nb_ts=1
	list_ts = range(ts_start,ts_end,nb_ts) 

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
	print((U_cine**2)*200.e-6*Sur/2*(rhol+rhov)/2)

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

	#Calcul de dx côté vapeur
	dxv=0

	dxv=(x1-x0)
	dxvt=50*dt
	ddxvt=(dxv-dxvt)/dxvt*100

	#Vérification de dx côté liquide (inutile)
	dxl=0
	dxl=dxv
	dxlt=50*dt
	ddxlt=(dxl-dxlt)/dxlt*100

	#U est constant dans la phase liquude
	U=U_sor

	#Definition d'un volume de contrôle
	def xl(t,U):
		a=0
		a=U*t+10.80e-3
		return a


	#Energy rise in the liquid 
	DELT=0

	DELT=(Tl1-Tl0)*cpl*rhol
	ETL1=Tl1*rhol*cpl
	ETL0=Tl0*rhol*cpl
	Pui2=DELT/dt

	#Enrgy rise in the vapor
	DEVT=0

	DEVT=(Tv1-Tv0)*cpv*rhov
	ETV1=Tv1*rhov*cpv
	ETV0=Tv0*rhov*cpv
	Pui3=DEVT/dt

	#Latent Energy Liquid Side

	DELL=0

	DELL=mpl1*dt*Lv_l
	MpointLT=kl*gradTl/Lv_l
	DELLT=-MpointLT*Lv_l*dt*Sur
	DDELL=(DELL-DELLT)/DELLT*100
	Pui4=DELL/dt


	#Latent Energy Vapor Side

	DELV=0

	DELV=mpv1*dt*Lv_l
	MpointVT=kv*gradTv/Lv_l
	DELVT=MpointVT*Lv_l*dt*Sur
	#DDELV=(DELV-DELVT)/DELVT*100
	Pui5=DELV/dt
	#Latent Energy Vapor Side

	#DEC=conv_sortie*dt
	#DECT=inter*(1/rhol-1/rhov)*cpl*rhol*Tls1
	#DDECT=(DEC-DECT)/DECT*100

	#Enrgie Convection
	Ener_conv=U_sor*Sur*Tvs1*cpv*dt*rhov-U_ent*Sur*Tle1*cpl*dt*rhov
	Pui6=Ener_conv/dt

	'''#Enrgie Gradient Pression
	Pre_Ener=PP*dt
	Pui7=PP'''

	#Calcul de la vitesse de l'interface
	Vit_i=(x1-x0)/dt

	#Défintion de x(t)
	def x(t,Vit_i):
		a=0
		a=Vit_i*t+10.40e-3
		return a

	#Bilan de puissances côté vapeur
	Pui_vap=-Pui1+Pui3-Pui5#+Pui12

	#Bilan de puissances côté liquide
	Pui_liq=Pui2-Pui4+Pui6

	#Bilan d'énergie
	Ener_Sum=DEVT+DELT-DELV-DELL+Ener_conv-diff_sortie*dt
	#-Pre_Ener
	Pui_Sum=Ener_Sum/dt





	#f.write("# Comments  \n")
	#f.write("# Legend\t t0=0.\t\t t1=1.e-8\t  Delta\t\t Err \n")

	f = open("Interface_Position_v_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxv,ddxvt)
	f.write(st)
	f.close()
	f = open("Interface_Position_l_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxl,ddxlt)
	f.write(st)
	f.close()
	f = open("T_Ener_Vap_vitesse.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETV0,ETV1,DEVT)
	f.write(st)
	f.close()
	f = open("T_Ener_Liq_vitesse.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETL0,ETL1,DELT)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Vap_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELV,DELV,100)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Liq_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELL,DELL,DDELL)
	f.write(st)
	f.close()
	f = open("Ener_Conv_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,Ener_conv,Ener_conv)
	f.write(st)
	f.close()
	f = open("Flux_limite_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,diff_sortie*dt,diff_sortie*dt)
	f.write(st)
	f.close()
	'''f = open("Pressure_Gradient_vitesse.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,Pre_Ener,Pre_Ener)
	f.write(st)
	f.close()'''
	f = open("Delta_Energy_Sum_vitesse.txt", "w")
	st= "-\t\t -\t %g\t - \n" %(Ener_Sum)
	f.write(st)
	f.close()
	f = open("Power_Sum_vitesse.txt", "w")
	st= "-\t\t -\t %g\t - \n" %(Pui_Sum)
	f.write(st)
	f.close()
	
	DeleteActivePlots()
	
	CloseDatabase("adiab/lata/post.lata")
	OpenDatabase("same_rho/lata/post.lata")

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
	ChangeActivePlotsVar("VITESSE_X_ELEM_dom")
	Query("Average value")
	U_cine=GetQueryOutputValue()

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

	f = open("same_rho/evevap_remesh_same_rho_pb_Diffusion_chaleur.out")
	lines=f.readlines()
	f.close()
	st = lines[-1]
	l=st.split()
	tps=float(l[0])
	diff_entree=float(l[1])
	Pui1=diff_entree

	inter=(kv*gradTv-kl*gradTl)/Lv_l

	#Calcul de dx côté vapeur
	dxv=0

	dxv=(x1-x0)
	dxvt=-(inter)/rhov*dt
	ddxvt=(dxv-dxvt)/dxvt*100

	#Vérification de dx côté liquide (inutile)
	dxl=0
	dxl=dxv
	dxlt=-(inter)/rhol*dt+inter*(1/rhol-1/rhov)*dt
	ddxlt=(dxl-dxlt)/dxlt*100

	#Definition d'un volume de contrôle
	def xl(t,U):
		a=0
		a=U*t+10.80e-3
		return a


	#Energy rise in the liquid 
	DELT=0

	DELT=(Tl1-Tl0)*cpl*rhol
	ETL1=Tl1*rhol*cpl
	ETL0=Tl0*rhol*cpl
	Pui2=DELT/dt

	#Enrgy rise in the vapor
	DEVT=0

	DEVT=(Tv1-Tv0)*cpv*rhov
	ETV1=Tv1*rhov*cpv
	ETV0=Tv0*rhov*cpv
	Pui3=DEVT/dt

	#Latent Energy Liquid Side

	DELL=0

	DELL=mpl1*dt*Lv_l
	MpointLT=kl*gradTl/Lv_l
	DELLT=-MpointLT*Lv_l*dt*Sur
	DDELL=(DELL-DELLT)/DELLT*100
	Pui4=DELL/dt

	#Latent Energy Vapor Side

	DELV=0

	DELV=mpv1*dt*Lv_l
	MpointVT=kv*gradTv/Lv_l
	DELVT=MpointVT*Lv_l*dt*Sur
	#DDELV=(DELV-DELVT)/DELVT*100
	Pui5=DELV/dt

	#Latent Energy Vapor Side

	#DEC=conv_sortie*dt
	#DECT=inter*(1/rhol-1/rhov)*cpl*rhol*Tls1
	#DDECT=(DEC-DECT)/DECT*100

	#Enrgie Convection
	Ener_conv=U_sor*Sur*Tvs1*cpv*dt*rhov-U_ent*Sur*Tle1*cpl*dt*rhol
	Pui6=Ener_conv/dt
	
	#Calcul de la vitesse de l'interface
	Vit_i=(x1-x0)/dt

	#Défintion de x(t)
	def x(t,Vit_i):
		a=0
		a=Vit_i*t+10.40e-3
		return a

	#Bilan de puissances côté vapeur
	Pui_vap=-Pui1+Pui3-Pui5

	#Bilan de puissances côté liquide
	Pui_vap=Pui2-Pui4 

	#Bilan d'énergie
	Ener_Sum=DEVT+DELT-DELV-DELL-diff_entree*dt+Ener_conv
	#-Pre_Ener
	Pui_Sum=Ener_Sum/dt

	#Bilan_phase_vapeur
	#Ener_Vap=DEVT+DELV-diff_entree/2.*dt



	#f.write("# Comments  \n")
	#f.write("# Legend\t t0=0.\t\t t1=1.e-8\t  Delta\t\t Err \n")

	f = open("Interface_Position_v_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxv,ddxvt)
	f.write(st)
	f.close()
	f = open("Interface_Position_l_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxl,ddxlt)
	f.write(st)
	f.close()
	f = open("T_Ener_Vap_same_rho.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETV0,ETV1,DEVT)
	f.write(st)
	f.close()
	f = open("T_Ener_Liq_same_rho.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETL0,ETL1,DELT)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Vap_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELV,DELV,0.)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Liq_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELL,DELL,DDELL)
	f.write(st)
	f.close()
	f = open("Flux_limite_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,diff_entree*dt,diff_entree*dt)
	f.write(st)
	f.close()
	f = open("Ener_Conv_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,Ener_conv,Ener_conv)
	f.write(st)
	f.close()
	'''f = open("Bilan_Ener_Vapeur_same_rho.txt", "w")
	st= "-\t\t -\t %g\t - \n" %(Ener_Vap)
	f.write(st)
	f.close()'''
	'''f = open("Pressure_Gradient_same_rho.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,Pre_Ener,Pre_Ener)
	f.write(st)
	f.close()'''
	f = open("Delta_Energy_Sum_same_rho.txt", "w")
	st= "-\t\t -\t %g\t - \n" %(Ener_Sum)
	f.write(st)
	f.close()
	f = open("Power_Sum_same_rho.txt", "w")
	st= "-\t-\t%g\t-\n" %(Pui_Sum)
	f.write(st)
	f.close()

	DeleteActivePlots()

	CloseDatabase("same_rho/lata/post.lata")
	OpenDatabase("std_2pas/lata/post.lata")
	
	dt=dt2
	rhov=rhol/2
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
	SetTimeSliderState(1)
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
	ChangeActivePlotsVar("VITESSE_X_ELEM_dom")
	Query("Average value")
	U_cine=GetQueryOutputValue()


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
	
	f = open("std_2pas/evap_remesh_2_pas_pb_Diffusion_chaleur.out")
	lines=f.readlines()
	f.close()
	st = lines[-1]
	l=st.split()
	tps=float(l[0])
	diff_entree=float(l[1])
	Pui1=diff_entree

	inter=(kv*gradTv-kl*gradTl)/Lv_l

	#Calcul de dx côté vapeur
	dxv=0

	dxv=(x1-x0)
	dxvt=-(inter)/rhov*dt
	ddxvt=(dxv-dxvt)/dxvt*100
	print dxvt
	
	#Vérification de dx côté liquide (inutile)
	dxl=0
	dxl=dxv
	dxlt=-(inter)/rhol*dt+inter*(1/rhol-1/rhov)*dt
	ddxlt=(dxl-dxlt)/dxlt*100
	print dxlt

	#Definition d'un volume de contrôle
	def xl(t,U):
		a=0
		a=U*t+10.80e-3
		return a


	#Energy rise in the liquid 
	DELT=0

	DELT=(Tl1-Tl0)*cpl*rhol
	ETL1=Tl1*rhol*cpl
	ETL0=Tl0*rhol*cpl
	Pui2=DELT/dt

	#Enrgy rise in the vapor
	DEVT=0

	DEVT=(Tv1-Tv0)*cpv*rhov
	ETV1=Tv1*rhov*cpv
	ETV0=Tv0*rhov*cpv
	Pui3=DEVT/dt

	#Latent Energy Liquid Side

	DELL=0

	DELL=mpl1*dt*Lv_l
	MpointLT=kl*gradTl/Lv_l
	DELLT=-MpointLT*Lv_l*dt*Sur
	DDELL=(DELL-DELLT)/DELLT*100
	Pui4=DELL/dt

	#Latent Energy Vapor Side

	DELV=0

	DELV=mpv1*dt*Lv_l
	MpointVT=kv*gradTv/Lv_l
	DELVT=MpointVT*Lv_l*dt*Sur
	#DDELV=(DELV-DELVT)/DELVT*100
	Pui5=DELV/dt

	#Latent Energy Vapor Side

	#DEC=conv_sortie*dt
	#DECT=inter*(1/rhol-1/rhov)*cpl*rhol*Tls1
	#DDECT=(DEC-DECT)/DECT*100

	#Calcul de la vitesse de l'interface
	Vit_i=(x1-x0)/dt

	#Défintion de x(t)
	def x(t,Vit_i):
		a=0
		a=Vit_i*t+10.40e-3
		return a

	#Enrgie Convection
	Ener_conv=U_sor*Sur*Tvs1*cpv*dt*rhov-U_ent*Sur*Tle1*cpl*dt*rhol
	Pui6=Ener_conv/dt
	
	#Bilan de puissances côté vapeur
	Pui_vap=-Pui1+Pui3-Pui5

	#Bilan de puissances côté liquide
	Pui_vap=Pui2-Pui4 

	#Bilan d'énergie
	Ener_Sum=DEVT+DELT-DELV-DELL-diff_entree*dt+Ener_conv
	#-Pre_Ener
	Pui_Sum=Ener_Sum/dt

	#Bilan_phase_vapeur
	#Ener_Vap=DEVT+DELV-diff_entree/2.*dt



	#f.write("# Comments  \n")
	#f.write("# Legend\t t0=0.\t\t t1=1.e-8\t  Delta\t\t Err \n")

	f = open("Interface_Position_v_2_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxv,ddxvt)
	f.write(st)
	f.close()
	f = open("Interface_Position_l_2_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxl,ddxlt)
	f.write(st)
	f.close()
	f = open("T_Ener_Vap_2_pas.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETV0,ETV1,DEVT)
	f.write(st)
	f.close()
	f = open("T_Ener_Liq_2_pas.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETL0,ETL1,DELT)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Vap_2_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELV,DELV,0.)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Liq_2_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELL,DELL,DDELL)
	f.write(st)
	f.close()
	f = open("Flux_limite_2_pas.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,diff_entree*dt,diff_entree*dt)
	f.write(st)
	f.close()
	f = open("Ener_Conv_2_pas.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,Ener_conv,Ener_conv)
	f.write(st)
	f.close()
	f = open("Delta_Energy_Sum_2_pas.txt", "w")
	st= "-\t\t -\t %g\t - \n" %(Ener_Sum)
	f.write(st)
	f.close()
	f = open("Power_Sum_2_pas.txt", "w")
	st= "-\t-\t%g\t-\n" %(Pui_Sum)
	f.write(st)
	f.close()
	
	dt=dt_sauv
	
	DeleteActivePlots()
	
	CloseDatabase("std_2pas/lata/post.lata")
	OpenDatabase("std_30pas/lata/post.lata")
	
	dt=30*dt/9
	rhov=rhol/2
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
	SetTimeSliderState(30)
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
	ChangeActivePlotsVar("VITESSE_X_ELEM_dom")
	Query("Average value")
	U_cine=GetQueryOutputValue()


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
	
	f = open("std_30pas/evap_remesh_30_pas_pb_Diffusion_chaleur.out")
	lines=f.readlines()
	f.close()
	st = lines[-1]
	l=st.split()
	tps=float(l[0])
	diff_entree=float(l[1])
	Pui1=diff_entree

	inter=(kv*gradTv-kl*gradTl)/Lv_l

	#Calcul de dx côté vapeur
	dxv=0

	dxv=(x1-x0)
	dxvt=-(inter)/rhov*dt
	ddxvt=(dxv-dxvt)/dxvt*100
	print dxvt
	
	#Vérification de dx côté liquide (inutile)
	dxl=0
	dxl=dxv
	dxlt=-(inter)/rhol*dt+inter*(1/rhol-1/rhov)*dt
	ddxlt=(dxl-dxlt)/dxlt*100
	print dxlt

	#Definition d'un volume de contrôle
	def xl(t,U):
		a=0
		a=U*t+10.80e-3
		return a


	#Energy rise in the liquid 
	DELT=0

	DELT=(Tl1-Tl0)*cpl*rhol
	ETL1=Tl1*rhol*cpl
	ETL0=Tl0*rhol*cpl
	Pui2=DELT/dt

	#Enrgy rise in the vapor
	DEVT=0

	DEVT=(Tv1-Tv0)*cpv*rhov
	ETV1=Tv1*rhov*cpv
	ETV0=Tv0*rhov*cpv
	Pui3=DEVT/dt

	#Latent Energy Liquid Side

	DELL=0

	DELL=mpl1*dt*Lv_l
	MpointLT=kl*gradTl/Lv_l
	DELLT=-MpointLT*Lv_l*dt*Sur
	DDELL=(DELL-DELLT)/DELLT*100
	Pui4=DELL/dt

	#Latent Energy Vapor Side

	DELV=0

	DELV=mpv1*dt*Lv_l
	MpointVT=kv*gradTv/Lv_l
	DELVT=MpointVT*Lv_l*dt*Sur
	#DDELV=(DELV-DELVT)/DELVT*100
	Pui5=DELV/dt

	#Latent Energy Vapor Side

	#DEC=conv_sortie*dt
	#DECT=inter*(1/rhol-1/rhov)*cpl*rhol*Tls1
	#DDECT=(DEC-DECT)/DECT*100

	#Calcul de la vitesse de l'interface
	Vit_i=(x1-x0)/dt

	#Défintion de x(t)
	def x(t,Vit_i):
		a=0
		a=Vit_i*t+10.40e-3
		return a

	#Enrgie Convection
	Ener_conv=U_sor*Sur*Tvs1*cpv*dt*rhov-U_ent*Sur*Tle1*cpl*dt*rhol
	Pui6=Ener_conv/dt
	
	#Bilan de puissances côté vapeur
	Pui_vap=-Pui1+Pui3-Pui5

	#Bilan de puissances côté liquide
	Pui_vap=Pui2-Pui4 

	#Bilan d'énergie
	Ener_Sum=DEVT+DELT-DELV-DELL-diff_entree*dt+Ener_conv
	#-Pre_Ener
	Pui_Sum=Ener_Sum/dt

	#Bilan_phase_vapeur
	#Ener_Vap=DEVT+DELV-diff_entree/2.*dt



	#f.write("# Comments  \n")
	#f.write("# Legend\t t0=0.\t\t t1=1.e-8\t  Delta\t\t Err \n")

	f = open("Interface_Position_v_30_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxv,ddxvt)
	f.write(st)
	f.close()
	f = open("Interface_Position_l_30_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(x0,x1,dxl,ddxlt)
	f.write(st)
	f.close()
	f = open("T_Ener_Vap_30_pas.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETV0,ETV1,DEVT)
	f.write(st)
	f.close()
	f = open("T_Ener_Liq_30_pas.txt", "w")
	st= "%g\t %g\t %g\t - \n" %(ETL0,ETL1,DELT)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Vap_30_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELV,DELV,0.)
	f.write(st)
	f.close()
	f = open("Lat_Ener_Liq_30_pas.txt", "w")
	st= "%g\t\t %g\t %g\t %g \n" %(0.0,DELL,DELL,DDELL)
	f.write(st)
	f.close()
	f = open("Flux_limite_30_pas.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,diff_entree*dt,diff_entree*dt)
	f.write(st)
	f.close()
	f = open("Ener_Conv_30_pas.txt", "w")
	st= "%g\t\t %g\t %g\t - \n" %(0.0,Ener_conv,Ener_conv)
	f.write(st)
	f.close()
	f = open("Delta_Energy_Sum_30_pas.txt", "w")
	st= "-\t\t -\t %g\t - \n" %(Ener_Sum)
	f.write(st)
	f.close()
	f = open("Power_Sum_30_pas.txt", "w")
	st= "-\t-\t%g\t-\n" %(Pui_Sum)
	f.write(st)
	f.close()
	
	quit()
