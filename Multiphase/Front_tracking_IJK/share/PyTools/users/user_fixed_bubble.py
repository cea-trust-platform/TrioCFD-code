# -*- coding: utf-8 -*-
#import rlcompleter
from Tools_fixedBubbles import *
import readline
import unittest
import pickle
import os


from matplotlib.pyplot import clf
#from scipy.sparse.linalg.isolve.minres import Ainfo

readline.parse_and_bind('tab:complete')

#from cmath import sin
#from scipy import misc
#from matplotlib.axes._subplots import Subplot
# import Tools
#import numpy as np
#import matplotlib.pyplot as plt
clf
global compteur
compteur=0
os.system('rm -f *.png')


def Init_case() : 
    res = {}
    return res

def Add_section(dic, name) :
    section = {}
    dic[name] = section
    return 

def UpdateParametricVariables(Run) :
    

    

        if (Run["variable"]["fstatBIT"]!=None):
		Run["parametric"]["Reb"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["mu"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 11 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 22 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 33 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 11 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 22 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 33 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["production 11 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["production 11 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 11 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 22 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 33 WIF"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 11 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 22 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 33 WIT"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 11 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 22 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 33 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 11 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 22 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 33 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["production 11 tot fixe"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["diffusion_3_1"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["diffusion_3_2"][0].append(Run['parametric']["parametre"])
	else:
		Run["parametric"]["epsilon 11 tot libre"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 22 tot libre"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["epsilon 33 tot libre"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 11 tot libre"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 22 tot libre"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["redistribution 33 tot libre"][0].append(Run['parametric']["parametre"])
		Run["parametric"]["production 11 tot libre"][0].append(Run['parametric']["parametre"])


        if (Run["variable"]["fstatBIT"]!=None):
		fstat=Run["variable"]["fstatBIT"]
		i=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
		iw_int=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IUint"],r'$U_l$', r'$z$', r'$u_l$')
		iww_int=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTUI_INT"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
		iw_intww=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTUIUI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$') 	
                iwww_int=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTUI_INTUI_INT"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')

		fstat=Run["variable"]["fstat"]
		iw=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI"],r'$U_l$', r'$z$', r'$u_l$')
		iww=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
		iwww=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')  
	
               

		Run["parametric"]["Reb"][1].append(Run['parametric']["parametre"])
		Run["parametric"]["mu"][1].append(Run["variable"]["mu"])
		Run["parametric"]["epsilon 11 WIF"][1].append(Run['Rij_BIT']["dissipation"].compo(0,0).mean())
		Run["parametric"]["epsilon 22 WIF"][1].append(Run['Rij_BIT']["dissipation"].compo(1,1).mean())
		Run["parametric"]["epsilon 33 WIF"][1].append(Run['Rij_BIT']["dissipation"].compo(2,2).mean())
		Run["parametric"]["epsilon 11 WIT"][1].append(Run['Rij_SIT']["dissipation"].compo(0,0).mean())
		Run["parametric"]["epsilon 22 WIT"][1].append(Run['Rij_SIT']["dissipation"].compo(1,1).mean())
		Run["parametric"]["epsilon 33 WIT"][1].append(Run['Rij_SIT']["dissipation"].compo(2,2).mean())
		Run["parametric"]["redistribution 11 WIF"][1].append(Run['Rij_BIT']["redistribution"].compo(0,0).mean())
		Run["parametric"]["redistribution 22 WIF"][1].append(Run['Rij_BIT']["redistribution"].compo(1,1).mean())
		Run["parametric"]["redistribution 33 WIF"][1].append(Run['Rij_BIT']["redistribution"].compo(2,2).mean())
		Run["parametric"]["redistribution 11 WIT"][1].append(Run['Rij_SIT']["redistribution"].compo(0,0).mean())
		Run["parametric"]["redistribution 22 WIT"][1].append(Run['Rij_SIT']["redistribution"].compo(1,1).mean())
		Run["parametric"]["redistribution 33 WIT"][1].append(Run['Rij_SIT']["redistribution"].compo(2,2).mean())
		Run["parametric"]["production 11 WIF"][1].append(Run['Rij_BIT']["interface"].compo(0,0).mean())
		Run["parametric"]["production 11 WIT"][1].append(Run['Rij_SIT']["interface"].compo(0,0).mean())
		Run["parametric"]["epsilon 11 tot fixe"][1].append(Run['Rij_tot']["dissipation"].compo(0,0).mean())
		Run["parametric"]["epsilon 22 tot fixe"][1].append(Run['Rij_tot']["dissipation"].compo(1,1).mean())
		Run["parametric"]["epsilon 33 tot fixe"][1].append(Run['Rij_tot']["dissipation"].compo(2,2).mean())
		Run["parametric"]["redistribution 11 tot fixe"][1].append(Run['Rij_tot']["redistribution"].compo(0,0).mean())
		Run["parametric"]["redistribution 22 tot fixe"][1].append(Run['Rij_tot']["redistribution"].compo(1,1).mean())
		Run["parametric"]["redistribution 33 tot fixe"][1].append(Run['Rij_tot']["redistribution"].compo(2,2).mean())
		Run["parametric"]["production 11 tot fixe"][1].append(Run['Rij_tot']["interface"].compo(0,0).mean())
		Run["parametric"]["diffusion_3_1"][1].append(iwww.mean()*(-1.))
		Run["parametric"]["diffusion_3_2"][1].append(((iw_intww-iwww_int)*3.).mean()*(-1.))
	else:
		Run["parametric"]["epsilon 11 tot libre"][1].append(Run['Rij_tot']["dissipation"].compo(0,0).mean())
		Run["parametric"]["epsilon 22 tot libre"][1].append(Run['Rij_tot']["dissipation"].compo(1,1).mean())
		Run["parametric"]["epsilon 33 tot libre"][1].append(Run['Rij_tot']["dissipation"].compo(2,2).mean())
		Run["parametric"]["redistribution 11 tot libre"][1].append(Run['Rij_tot']["redistribution"].compo(0,0).mean())
		Run["parametric"]["redistribution 22 tot libre"][1].append(Run['Rij_tot']["redistribution"].compo(1,1).mean())
		Run["parametric"]["redistribution 33 tot libre"][1].append(Run['Rij_tot']["redistribution"].compo(2,2).mean())
		Run["parametric"]["production 11 tot libre"][1].append(Run['Rij_tot']["interface"].compo(0,0).mean())
	return

def CalculerReb(Run):
    fstat=Run["variable"]["fstat"]
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Iv=(I-1.0)*(-1.0)
    u_l=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')*I.inv()).compo(0)
    u_v=(Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_l$', r'$z$', r'$u_l$')*Iv.inv()).compo(0)
    Run["parametric"]["parametre"]=((u_v-u_l)*(Run["variable"]["db"]*Run["variable"]["rho_l"]/Run["variable"]["mu"])).mean()
    Run["parametric"]["cd"][1].append((4./3.)*Run["variable"]["alpha"]*Run["variable"]["db"]*(Run["variable"]["rho_l"]-Run["variable"]["rho_v"])*abs(Run["variable"]["gravity"])*(1./(((u_v-u_l).mean())**2)))
    Run["parametric"]["cd"][0].append(Run["parametric"]["parametre"])
    return




def UsualVariable_fixedBubble(Run) :
    fstat=Run["variable"]["fstat"]

    ################################################################################
    ########### PLOT DES PROFIL DE VITESSE ET INDICATRICE ##########################
    ################################################################################
    #sym_tens=[1,-1,1,-1,1,-1,1,-1,1]
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    xlims=[0,0.0015]


    
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_l$', r'$z$', r'$u_l$')
    u=u_l
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    



    if (Run["variable"]["fstatBIT"]!=None):
	    u_np=Field.LoadFromFile(fstat,["coordonnee_K"] , ["U_NOPERTURBE", "V_NOPERTURBE", "W_NOPERTURBE"],r'$U_l$', r'$z$', r'$u_l$')
	    u_l_int=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IUint", "IVint", "IWint"],r'$U_l$', r'$z$', r'$u_l$')
	    uu_l_spa=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTUI_INT", "UI_INTVI_INT", "UI_INTWI_INT","UI_INTVI_INT", "VI_INTVI_INT", "VI_INTWI_INT","UI_INTWI_INT", "VI_INTWI_INT", "WI_INTWI_INT"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')

    

    uu=uu_l
    ufu=u.MatProduct(u)
    Rij=uu-ufu

    if (Run["variable"]["fstatBIT"]!=None):
	    ufu_int=u_l_int.MatProduct(u_l_int)
	    Rij_spa=uu_l_spa-ufu_int
	    Rij_tmp=Rij-Rij_spa


    # Set names : 
    u_l.settex(r'$u_{liquide}$ spherical (TrioCFD)')
    Iv=(I-1.0)*(-1.0)


    ############################################################################################
    ############################## PLOT DES U_RMS ##############################################
    ############################################################################################

    #tracer(Run,[Rij.compo(0,0), Rij_spa.compo(0,0), Rij_tmp.compo(0,0)], 'uu',textitle=r"Axial $\overline{u'u'}$", couleur=[3,1])

    #tracer(Run,[Rij.compo(2,2), Rij_spa.compo(2,2), Rij_tmp.compo(2,2)], 'vv',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[3,1])

    #tracer(Run,[Rij.compo(1,1), Rij_spa.compo(1,1), Rij_tmp.compo(1,1)], 'ww',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[3,1])
    u_l=u_l*I.inv()
    u_v=u_v*Iv.inv()
    Run["variable"]['I']=Iv
    Run["variable"]['u_l']=u_l
    Run["variable"]['u_v']=u_v

    if (Run["variable"]["fstatBIT"]!=None):
    	Run["variable"]['Rij_wit']=Rij_tmp
    if (Run["variable"]["fstatBIT"]!=None):
	Run["variable"]['Rij_wif']=Rij_spa

    if (Run["option"]["debug"]) :
		tracer(Run,[u_l, u_v], 'u_l_v',  [0])


    if (Run["variable"]["fstatBIT"]!=None):
	    Run["parametric"]["R11 WIF"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R22 WIF"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R33 WIF"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R11 WIT"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R22 WIT"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R33 WIT"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R11 WIF"][1].append((Rij_spa.compo(0,0)._npa[0,:]).mean())
	    Run["parametric"]["R22 WIF"][1].append((Rij_spa.compo(1,1)._npa[0,:]).mean())
	    Run["parametric"]["R33 WIF"][1].append((Rij_spa.compo(2,2)._npa[0,:]).mean())
	    Run["parametric"]["R11 WIT"][1].append((Rij_tmp.compo(0,0)._npa[0,:]).mean())
	    Run["parametric"]["R22 WIT"][1].append((Rij_tmp.compo(1,1)._npa[0,:]).mean())
	    Run["parametric"]["R33 WIT"][1].append((Rij_tmp.compo(2,2)._npa[0,:]).mean())
	    Run["parametric"]["uvfixe"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["ulfixe"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["urfixe"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["uvfixe"][1].append(((u_v).compo(0)._npa[0,:]).mean())
	    Run["parametric"]["ulfixe"][1].append(((u_l).compo(0)._npa[0,:]).mean())
	    Run["parametric"]["urfixe"][1].append(((u_v-u_l).compo(0)._npa[0,:]).mean())

    else : 
	    Run["parametric"]["R11 tot"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R22 tot"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R33 tot"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["R11 tot"][1].append((Rij.compo(0,0)._npa[0,:]).mean())
	    Run["parametric"]["R22 tot"][1].append((Rij.compo(1,1)._npa[0,:]).mean())
	    Run["parametric"]["R33 tot"][1].append((Rij.compo(2,2)._npa[0,:]).mean())
	    Run["parametric"]["uvlibre"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["ullibre"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["urlibre"][0].append(Run["parametric"]["parametre"])
	    Run["parametric"]["uvlibre"][1].append(((u_v).compo(0)._npa[0,:]).mean())
	    Run["parametric"]["ullibre"][1].append(((u_l).compo(0)._npa[0,:]).mean())
	    Run["parametric"]["urlibre"][1].append(((u_v-u_l).compo(0)._npa[0,:]).mean())




    #print "le rapport BIT/SIT en xx est ", (Rij_spa.compo(0,0)._npa[0,:]).mean()/(Rij_tmp.compo(0,0)._npa[0,:]).mean()
    #print "le rapport BIT/SIT en yy est ", (Rij_spa.compo(1,1)._npa[0,:]).mean()/(Rij_tmp.compo(1,1)._npa[0,:]).mean()
    #print "le rapport BIT/SIT en zz est ", (Rij_spa.compo(2,2)._npa[0,:]).mean()/(Rij_tmp.compo(2,2)._npa[0,:]).mean()
    #print "le rapport BIT/SIT en xz est ", (Rij_spa.compo(2,0)._npa[0,:]).mean()/(Rij_tmp.compo(2,0)._npa[0,:]).mean()
    return

############################################################################################
############################## CHARGEMENT DES CHAMPS #######################################
############################################################################################



def RunDiphasicCase_fixedBubble(Run):
	gravity=Run["variable"]["gravity"]
	rho_l=Run["variable"]["rho_l"]
	rho_v=Run["variable"]["rho_v"]
	alpha_v=Run["variable"]["alpha"]
	mu=Run["variable"]["mu"]
	sigma=Run["variable"]["sigma"]
	Source=(Run["variable"]["retau"]*Run["variable"]["mu"])**2
	fstat=Run["variable"]["fstat"]
	uwall=Run["variable"]["u_inf"]
	dvar = Field.getEntries(fstat)
    #################################### DONNEES RIJ TOTAL 
	if (Run["option"]["Rij_tot"]) :
		fstat=Run["variable"]["fstat"] 
		I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')

		if (Run["variable"]["fstatBIT"]!=None) :
			I_np=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I_NP"],'I', 'z', r'I')

		uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
		uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
		u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
		dudxdudx_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$', bc=2)
		dudydudy_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$', bc=2)
		dudzdudz_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$', bc=2)	    

		## les termes avec des pressions
		if (Run["option"]["pression_NP_pour_tot"] and Run["variable"]["fstatBIT"]!=None) :
			pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_NOPERTURBE"],'P_l', 'z', r'$P_l$')*(I_np.inv(0.000001))*I
			pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIUNP", "PIVNP", "PIWNP"],r'$Up_l$', r'$z$', r'$up_l$')*(I_np.inv(0.000001))*I
			pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdxNP","IPdUdyNP","IPdUdzNP","IPdVdxNP", "IPdVdyNP","IPdVdzNP","IPdWdxNP","IPdWdyNP","IPdWdzNP"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')*(I_np.inv(0.000001))*I
		else :
			# avec les pressions réelles
			pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
			pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
			pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')

		Rij_l_TOT=fluctuij(u_l, u_l, uu_l)
		pdu_ll_TOT=fluctuij(pression_l, u_l.grad(2), pdu_l)
		pu_ll_TOT=fluctuij(pression_l, u_l, pu_l)
		dudu0_TOT=dudxdudx_l+dudydudy_l+dudzdudz_l
		dudu_l_TOT=fluctuij(u_l.grad(2), u_l.grad(2), dudu0_TOT)
		Rijk_l_TOT=fluctuijk(uuu_l, u_l, uu_l)         
	  	u_l_TOT=u_l*1.0



    #################################### DONNEES RIJ BIT 
	if (Run["option"]["Rij_BIT"] and Run["variable"]["fstatBIT"]!=None) :
		fstat=Run["variable"]["fstatBIT"] 
		I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
		I_np=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I_NP"],'I', 'z', r'I')
		uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTUI_INT", "UI_INTVI_INT", "UI_INTWI_INT","UI_INTVI_INT", "VI_INTVI_INT", "VI_INTWI_INT","UI_INTWI_INT", "VI_INTWI_INT", "WI_INTWI_INT"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
		uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTUI_INTUI_INT", "UI_INTUI_INTVI_INT", "UI_INTUI_INTWI_INT","UI_INTUI_INTVI_INT", "UI_INTVI_INTVI_INT", "UI_INTVI_INTWI_INT","UI_INTUI_INTWI_INT", "UI_INTVI_INTWI_INT", "UI_INTWI_INTWI_INT","UI_INTUI_INTVI_INT", "UI_INTVI_INTVI_INT", "UI_INTVI_INTWI_INT","UI_INTVI_INTVI_INT", "VI_INTVI_INTVI_INT", "VI_INTVI_INTWI_INT","UI_INTVI_INTWI_INT", "VI_INTVI_INTWI_INT", "VI_INTWI_INTWI_INT","UI_INTUI_INTWI_INT", "UI_INTVI_INTWI_INT", "UI_INTWI_INTWI_INT","UI_INTVI_INTWI_INT", "VI_INTVI_INTWI_INT", "VI_INTWI_INTWI_INT","UI_INTWI_INTWI_INT", "VI_INTWI_INTWI_INT", "WI_INTWI_INTWI_INT" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
		u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IUint", "IVint", "IWint"],r'$U_l$', r'$z$', r'$u_l$')
		dudxdudx_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUINTdxdUINTdx", "dUINTdxdVINTdx", "dUINTdxdWINTdx","dUINTdxdVINTdx", "dVINTdxdVINTdx", "dVINTdxdWINTdx","dUINTdxdWINTdx", "dVINTdxdWINTdx", "dWINTdxdWINTdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$', bc=2)
		dudydudy_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUINTdydUINTdy", "dUINTdydVINTdy", "dUINTdydWINTdy","dUINTdydVINTdy", "dVINTdydVINTdy", "dVINTdydWINTdy","dUINTdydWINTdy", "dVINTdydWINTdy", "dWINTdydWINTdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$', bc=2)
		dudzdudz_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUINTdzdUINTdz", "dUINTdzdVINTdz", "dUINTdzdWINTdz","dUINTdzdVINTdz", "dVINTdzdVINTdz", "dVINTdzdWINTdz","dUINTdzdWINTdz", "dVINTdzdWINTdz", "dWINTdzdWINTdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$', bc=2)	   
		# avec les pressions réelles
		pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IP_INT"],'P_l', 'z', r'$P_l$') #### a modifier par le Pint !!!!!!
		pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI_INTP_INT", "VI_INTP_INT", "WI_INTP_INT"],r'$Up_l$', r'$z$', r'$up_l$')
		pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"], ["P_INTdUINTdx","P_INTdUINTdy","P_INTdUINTdz","P_INTdVINTdx","P_INTdVINTdy","P_INTdVINTdz","P_INTdWINTdx","P_INTdWINTdy","P_INTdWINTdz"], 'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')

		Rij_l_BIT=fluctuij(u_l, u_l, uu_l)
		pdu_ll_BIT=fluctuij(pression_l, u_l.grad(2), pdu_l)
		pu_ll_BIT=fluctuij(pression_l, u_l, pu_l)
		dudu0_BIT=dudxdudx_l+dudydudy_l+dudzdudz_l
		dudu_l_BIT=fluctuij(u_l.grad(2), u_l.grad(2), dudu0_BIT)
		Rijk_l_BIT=fluctuijk(uuu_l, u_l, uu_l)         
		u_l_BIT=u_l*1.0


    #################################### DONNEES RIJ SIT

	if (Run["option"]["Rij_SIT"] and Run["variable"]["fstatBIT"]!=None) :	
		Rij_l_SIT  =  Rij_l_TOT- Rij_l_BIT
		pdu_ll_SIT = pdu_ll_TOT-pdu_ll_BIT
		pu_ll_SIT  =  pu_ll_TOT- pu_ll_BIT
		dudu_l_SIT = dudu_l_TOT-dudu_l_BIT
		u_l_SIT    =    u_l_TOT-   u_l_BIT
		Rijk_l_SIT = Rijk_l_TOT-Rijk_l_BIT
		if (Run["option"]["pression_NP_pour_SIT"]) :
			u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
			pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_NOPERTURBE"],'P_l', 'z', r'$P_l$')*(I_np.inv(0.000001))*I
			pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIUNP", "PIVNP", "PIWNP"],r'$Up_l$', r'$z$', r'$up_l$')*(I_np.inv(0.000001))*I
			pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdxNP","IPdUdyNP","IPdUdzNP","IPdVdxNP", "IPdVdyNP","IPdVdzNP","IPdWdxNP","IPdWdyNP","IPdWdzNP"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')*(I_np.inv(0.000001))*I	
			### correction pression SIT avec chmps non perturbe
			pdu_ll_SIT=fluctuij(pression_l, u_l.grad(2), pdu_l)
			pu_ll_SIT=fluctuij(pression_l, u_l, pu_l)

			### correction pression BIT avec nouveau SIT
			pdu_ll_BIT = pdu_ll_TOT-pdu_ll_SIT
			pu_ll_BIT = pu_ll_TOT-pu_ll_SIT     

	#################################### CALCUL RIJ_TOT _BIT et _SIT   



	if (Run["option"]["debug"] and Run["variable"]["fstatBIT"]!=None) :
		tracer(Run,[u_l_TOT,u_l_BIT,u_l_SIT], 'debug_ulx',  [0])
		tracer(Run,[u_l_TOT,u_l_BIT,u_l_SIT], 'debug_uly',  [1])
		tracer(Run,[u_l_TOT,u_l_BIT,u_l_SIT], 'debug_ulz',  [2])
		tracer(Run,[Rij_l_TOT, Rij_l_BIT, Rij_l_SIT], 'debug_R11',  [0])
		tracer(Run,[Rij_l_TOT, Rij_l_BIT, Rij_l_SIT], 'debug_R22',  [4])
		tracer(Run,[Rij_l_TOT, Rij_l_BIT, Rij_l_SIT], 'debug_R33',  [8])
		tracer(Run,[Rij_l_TOT, Rij_l_BIT, Rij_l_SIT], 'debug_R13',  [2])
		tracer(Run,[Rijk_l_TOT, Rijk_l_BIT, Rijk_l_SIT], 'debug_R111',  [0])
		tracer(Run,[pu_ll_TOT,pu_ll_BIT,pu_ll_SIT], 'debug_pulx',  [0])
		tracer(Run,[pu_ll_TOT,pu_ll_BIT,pu_ll_SIT], 'debug_puly',  [1])
		tracer(Run,[pu_ll_TOT,pu_ll_BIT,pu_ll_SIT], 'debug_pulz',  [2])
		tracer(Run,[pdu_ll_TOT,pdu_ll_BIT,pdu_ll_SIT], 'debug_pdu11',  [0])
		tracer(Run,[pdu_ll_TOT,pdu_ll_BIT,pdu_ll_SIT], 'debug_pdu22',  [4])
		tracer(Run,[pdu_ll_TOT,pdu_ll_BIT,pdu_ll_SIT], 'debug_pdu33',  [8])
		tracer(Run,[pdu_ll_TOT,pdu_ll_BIT,pdu_ll_SIT], 'debug_pdu13',  [2])
		tracer(Run,[dudu_l_TOT, dudu_l_BIT, dudu_l_SIT], 'debug_dudu11',  [0])
		tracer(Run,[dudu_l_TOT, dudu_l_BIT, dudu_l_SIT], 'debug_dudu22',  [4])
		tracer(Run,[dudu_l_TOT, dudu_l_BIT, dudu_l_SIT], 'debug_dudu33',  [8])
		tracer(Run,[dudu_l_TOT, dudu_l_BIT, dudu_l_SIT], 'debug_dudu13',  [2])
	if (Run["option"]["Rij_tot"]) :
	    Rij_diphase("Rij_tot", Run, u_l_TOT, Rij_l_TOT, Rijk_l_TOT, pu_ll_TOT, pdu_ll_TOT, rho_l, dudu_l_TOT, mu)
	if (Run["option"]["Rij_BIT"] and Run["variable"]["fstatBIT"]!=None) :
	    Rij_diphase("Rij_BIT", Run, u_l_BIT, Rij_l_BIT, Rijk_l_BIT, pu_ll_BIT, pdu_ll_BIT, rho_l, dudu_l_BIT, mu)
	if (Run["option"]["Rij_SIT"] and Run["option"]["Rij_BIT"] and Run["option"]["Rij_tot"] and Run["variable"]["fstatBIT"]!=None) :
	    Rij_diphase("Rij_SIT", Run, u_l_SIT, Rij_l_SIT, Rijk_l_SIT, pu_ll_SIT, pdu_ll_SIT, rho_l, dudu_l_SIT, mu)



	return



def Calcul_qdm_fixedbubble(Run):
    fstat=Run['variable']['fstat']
    gravity=Run['variable']['gravity']
    alpha_v=Run['variable']['alpha']
    rho_l=Run['variable']['rho_l']
    rho_v=Run['variable']['rho_v']
    mu=Run['variable']['mu']
    sigma=Run['variable']['sigma']
    Source=(Run['variable']['retau']*mu)**2
    dvar = Field.getEntries(fstat)
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii=Field.LoadFromFile(fstat,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    PaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PaiNx", "PaiNy", "PaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    aiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    vdvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","UdVdxaiNx", "VdVdxaiNx", "VdWdxaiNx","UdWdxaiNx", "VdWdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","UdVdyaiNy", "VdVdyaiNy", "VdWdyaiNy","UdWdyaiNy", "VdWdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","UdVdzaiNz", "VdVdzaiNz", "VdWdzaiNz","UdWdzaiNz", "VdWdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    rappel=Field.LoadFromFile(fstat,["coordonnee_K"] , ["Ivrappelx", "Ivrappelx", "Ivrappelx"],r'$U_l$', r'$z$', r'$u_l$')


    ############################################################################################
    ############################## PARAMETRES  #################################################
    ############################################################################################

    

    tracer(Run,[I], 'I')
    g=Field.initgravity([gravity,0,0], 'g', I)
    Iv=(I-1.0)*(-1.0)
    alpha_l=1-alpha_v
    rho_ll=rho_l*alpha_l
    rho_vv=rho_v*alpha_v
    rho_av=rho_ll+rho_vv
    S=Field.initsource([Source,0,0], 'beta', I)


    ############################################################################################
    ############################## BILAN DE QDM ################################################
    ############################################################################################
    
    rho1=I*rho_l
    rho2=Iv*rho_v
    rho=rho1+rho2
    u=u_v+u_l
    tracer(Run,[u],'u0', [0])
    tracer(Run,[u],'u1', [1])
    tracer(Run,[u],'u2', [2])
	
    ull=u_l*(I.inv(0.00001))
    uvv=u_v*(Iv.inv(0.00001))
   
    uu=uu_v+uu_l
    pression=pression_l+pression_v         
    rhog_av=g*rho_av
    rhog_tmp=rho*g
    rhog=rhog_tmp-rhog_av
    ufu=u.MatProduct(u)
    Rij=uu-ufu
    ufu_l=u_l.MatProduct(u_l)
    Rij_l=uu_l-ufu_l
    ufu_v=u_v.MatProduct(u_v)
    Rij_v=uu_v-ufu_v
    tracer(Run,[pression_l], 'pression')

    qdm_diph=qdm_diphase('diph',u, mu, pression, Rij, S, rho, aii, rho_av, g, sigma, I, rappel)
    qdm_liq=qdm_diphase('liquide',u_l, mu, pression_l, Rij_l, S, rho_l, aii, rho_av, g, sigma, I, rappel)
    qdm_gaz=qdm_diphase('gaz',u_v, mu, pression_v, Rij_v, S, rho_v, aii, rho_av, g, sigma, Iv, rappel)

    return


def Rij_diphase(name, Run, u, Rij, Rijk, pu, pdu, rho, dudu0, mu):

    rho_l=Run["variable"]["rho_l"]
    ### advection / production ###
    grad_Rij=Rij.grad(2)
    grad_u=u.grad(2)
    grad_uRij=((u.tensorProduct(Rij)).grad(2)).contract('ij','kijk')
    Rij_grad_u=(Rij.tensorProduct(grad_u)).contract('ij','kijk')
    Rij_grad_u_transpo=Rij_grad_u.transpo()

    advection=grad_uRij*(-1.0)
    production=(Rij_grad_u+Rij_grad_u_transpo)*(-1.0)


    #### diffusion ########
    grad_Rijk=(Rijk.grad(2)).contract('ij','kijk')
    diffusion=grad_Rijk*(-1.0)
    if Run["option"]["deactive_diffusion"]:
	diffusion=diffusion*0.0
    
    ####### transpression / redistribution ###############

    grad_pu=pu.grad(2)
    grad_pu_transpo=grad_pu.transpo()

    transpression = (grad_pu+grad_pu_transpo)*(-1.0)*(1./rho_l)
    redistribution = (pdu.transpo()+pdu)*(1./rho_l)

    #### dissipationn #######

    dissipation=dudu0*mu*2.0*(-1.0)*(1./rho_l)
    print "calcul dissipation rho_l=", rho_l, "et mu=", mu

    ####### diffusion moleculaire #########
    
    
    grad_grad_Rij=(grad_Rij.grad(2)).contract('ij', 'ijkk')
    diffusion_mol=grad_grad_Rij*mu*(1./rho_l)
    if Run["option"]["deactive_diffusion"]:
	diffusion_mol=diffusion_mol*0.0

    
    #### reconstruction theorique de la redistribution
    
    if Run["option"]["redistribution_theorique"]:
	    redistribution = (dissipation+diffusion_mol+diffusion-advection-production)*(-1.0)
	    redistribution._npa[0,0,:]=(redistribution.compo(1,1)*(-1.0))._npa[0,:]-(redistribution.compo(2,2))._npa[0,:]
	    transpression = transpression*0.0


    ##### production interfacial (residu)
    
    interface=(transpression+redistribution+dissipation+diffusion_mol+diffusion-advection-production)*(-1.0)


    ##### latex legend

    advection.settex(r'advection $\frac{DR_{i,j}}{Dt}$')
    production.settex(r'production $\overline{v^,_lv^,_j}\frac{\partial\overline{v_i}}{\partial x_l}+\overline{v^,_lv^,_i}\frac{\partial\overline{v_j}}{\partial x_l}$')
    redistribution.settex(r'redistribution $\overline{p^,\left(\frac{\partial v^,_j}{\partial x_i}+\frac{\partial v^,_i}{\partial x_j}\right)}$')
    dissipation.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_i}{\partial x_l} \frac{\partial v^,_j}{\partial x_l}}$')
    diffusion.settex(r'turbulent diffusion $\frac{\partial\overline{v^,_lv^,_iv^,_j}}{\partial x_l}$')
    transpression.settex(r'pressure diffusion $\frac{\partial}{\partial x_l}(\overline{p^,v^,_j}\delta_{il}+\overline{p^,v^,_i}\delta_{jl})$')
    diffusion_mol.settex(r'molecular diffusion  $\nu\frac{\partial^2\overline{v^,_iv^,_j}}{\partial^2 x_l}$')
    interface.settex(r'interface $\langle\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}\rangle+\langle\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I\rangle-\langle p^,(v^,_jn_i\delta^I+v^,_in_j\delta^I)\rangle$ ')
    

    ###################### plot ##################################

    if (Run["plot"]["plot_Rij"]) :
	    listing=[advection, production, redistribution, dissipation, diffusion, transpression, diffusion_mol, interface]
	    xlimm=Run["plot"]["xlim"]
	    tracer(Run,listing, name+"_11", [0])
	    tracer(Run,listing, name+"_22", [4])
	    tracer(Run,listing, name+"_33", [8])

	    tenseur_source_11=(redistribution+interface).compo(0,0)
	    tenseur_source_22=(redistribution+interface).compo(1,1)
	    tenseur_source_33=(redistribution+interface).compo(2,2)
	    tracer(Run,[tenseur_source_11,tenseur_source_22,tenseur_source_33], name+'_tenseur_source')
	    tracer(Run,[redistribution.compo(0,0),redistribution.compo(1,1),redistribution.compo(2,2)], name+'_redistribution')
	    tracer(Run,[interface.compo(0,0),interface.compo(1,1),interface.compo(2,2)], name+'_interface')

    Add_section(Run, name)
    Run[name]["advection"]=advection
    Run[name]["production"]=production
    Run[name]["redistribution"]=redistribution
    Run[name]["dissipation"]=dissipation
    Run[name]["diffusion"]=diffusion
    Run[name]["transpression"]=transpression
    Run[name]["diffusion_mol"]=diffusion_mol
    Run[name]["interface"]=interface




    return

def PlotDecompoCorrelationRij(Run, name):
        
	Run['Rij_SIT'][name].settex(name+' SIT')
	Run['Rij_BIT'][name].settex(name+' BIT')
	Run['Rij_tot'][name].settex(name+' TOT')

	tab=[Run['Rij_tot'][name], Run['Rij_BIT'][name], Run['Rij_SIT'][name]]
	tracer(Run,tab, name+'_decompo_11', [0])
	tracer(Run,tab, name+'_decompo_22', [4])
	tracer(Run,tab, name+'_decompo_33', [4])
	tracer(Run,tab, name+'_decompo_13', [6])



def qdm_diphase(phase, u, mu, pression, Rij, S, rho, aii=None, rho_av=None, g=None, sigma=None, I=None, rappel=None):
    """ 
    return a dictionary composed by all the qdm profil (turbulence, inertia; buoyancy etc...)
    @param the phase. Could be 'diph', 'liq' or 'gaz'
    @param the velocity field (vector field)
    @param viscosity scalar
    @param pressure field (scalar field)
    @param Reynolds stress tensor (tensor field)
    @param source (scalar field) 
    @param volumic mass (scalar)
    """
    

    tau1=u.grad(2)
    tau2=tau1+tau1.transpo()
    tau3=tau2.grad(2)
    tau4=tau3.contract('i', 'ikk')


    nonlin=u.tensorProduct(u)
    gradp=pression.grad(2)
    divuu=nonlin.grad(2)
    divuu=divuu.contract('i', 'ikk')
    divRij=Rij.grad(2)
    divRij=divRij.contract('i', 'ikk')
  
    
        ### rho et mu ponderation
    if phase=='diph':    
        divuu=divuu*rho
        divRij=rho*divRij
        divT=tau4*mu
        mutau=tau2*mu
    else:
        divuu=divuu*rho
        divRij=divRij*rho
        divT=tau4*mu
        mutau=tau2*mu       
  
    
    
    rhog_av=g*rho_av
    if phase=='diph':
        rhog_tmp=rho*g ## diphasique rho est un profil
    else:
        rhog_tmp=g*rho 
        
    if phase=='diph':
        rhog=rhog_tmp-rhog_av
        ai=aii*sigma
    else:
        rhog=rhog_tmp-rhog_av
        rhog=I*rhog
        S=I*S
        ai=divT+divRij*(-1.0)+gradp*(-1.0)+S+divuu*(-1.0)+rhog

    if 0:
	    divT.symetriser([1,1,1])
	    divRij.symetriser([1,1,1])
	    gradp.symetriser([1,1,1])
	    rhog.symetriser([1,1,1])
	    ai.symetriser([1,1,1])
    gradp=gradp*(-1.0)
    divRij=divRij*(-1.0)
    divuu=divuu*(-1.0)

    if phase=='diph':
    	insta=divT+divRij+gradp+S+divuu+rhog+ai
    else :
        insta=divT*0.0
    
    divuu.settex(r'inertia $-\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right)$') 
    divRij.settex(r'turbulence $-\nabla .\left(\overline{\rho\mathbf{v^,}\otimes\mathbf{v^,}}\right)$')     
    divT.settex(r'viscous stress $\nabla .\left(\nu\left(\nabla\overline{\mathbf{v}}+\nabla\overline{\mathbf{v}}^{t}\right)\right)$')  
    gradp.settex(r'pressure $-\nabla\overline{P}$')
    S.settex(r'source $-\beta\mathbf{e_x}$')
    rhog.settex(r'buoyancy $(\rho-\langle\rho\rangle)\mathbf{g}$')
    rappel.settex(r'rappel')
    insta.settex(r'residu')
         
    if phase=='diph':
        ai.settex(r'residue $\sigma\kappa\mathbf{n}\delta^i$')
    else:
        ai.settex(r'residue $\langle\frac{p}{\rho}\nabla\chi-\nu\nabla\mathbf{v}\nabla\chi\rangle$')

    if 1 :
        tracer(Run,[divT, divRij, gradp, S,rhog, ai, insta, rappel], phase+'_qdm_i', [0])
        tracer(Run,[divT, divRij, gradp, S,rhog, ai, insta, rappel], phase+'_qdm_j', [1])
        tracer(Run,[divT, divRij, gradp, S,rhog, ai, insta, rappel], phase+'_qdm_k', [2])
    
    res={'inertie':divuu, 'turbulence':divRij, 'viscous':divT, 'pression':gradp, 'source':S, 'interface':ai, 'rhog':rhog }
    
    stress_divRij=divRij.integ()
    stress_divuu=divuu.integ()
    stress_divT=divT.integ()
    stress_gradp=gradp.integ()
    stress_rhog=rhog.integ()
    stress_ai=ai.integ()
    stress_S=S.integ()
    stress_S=stress_S*(-1.0)
    stress_insta=insta.integ()
   

    stress_divuu.settex(r'inertia $\int^y_0\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right) dy$') 
    stress_divRij.settex(r'turbulence $\int^y_0\nabla .\left(\overline{\rho\mathbf{v^,}\otimes\mathbf{v^,}}\right) dy$')     
    stress_divT.settex(r'viscous stress $\int^y_0\nabla .\left(\nu\left(\nabla\overline{\mathbf{v}}+\nabla\overline{\mathbf{v}}^{t}\right)\right) dy$')  
    stress_gradp.settex(r'pressure $\int^y_0\nabla\overline{P}dy$')
    stress_S.settex(r'source $\int^y_0\beta dy$')
    stress_rhog.settex(r'buoyancy $\int^y_0(\rho-\langle\rho\rangle)\mathbf{g}dy$')     
    stress_ai.settex(r'interface $\int^y_0\sigma\kappa\mathbf{n}\delta^idy$')
    stress_insta.settex(r'residu')

    if phase=='diph':
        stress_ai.settex(r'interface (residue) $\int^y_0\sigma\kappa\mathbf{n}\delta^idy$')
    else:
        stress_ai.settex(r'interface (residue) $\int^y_0\langle\frac{p}{\rho}\nabla\chi-\nu\nabla\mathbf{v}\nabla\chi\rangle dy$')
    

    listing=[stress_divT, stress_divRij, stress_gradp, stress_S,stress_rhog, stress_ai, stress_insta]

    if 1 :   
        Run['plot']['zerocentre']=True
        tracer(Run,listing, phase+'_stress_i', [0])
        Run['plot']['zerocentre']=False
        tracer(Run,listing, phase+'_stress_j', [1])   
        tracer(Run,listing, phase+'_stress_k', [2])  
        
         
    #### Bilan des pression
    db=0.3
    p_iner=u.compo(0)*u.compo(0)*rho*(1./3.)*(-1.0)
    p_turb=(Rij.compo(0,0)+Rij.compo(1,1)+Rij.compo(2,2))*rho*(1./3.)*(-1.0)
    p_pi=rhog.compo(0)
    p_sig=(ai.compo(0)+ai.compo(1)+ai.compo(2))*(1./3.)
    p_tau=(mutau.compo(0,0)+mutau.compo(1,1)+mutau.compo(2,2))*(1./3.)
    p=p_turb+p_sig+p_tau+p_pi
    preal=pression
    p.settex(r'$p_{tot}$')
    preal.settex(r'$p_{DNS}$')
    p_iner.settex(r'$p_{inertie}$')
    p_turb.settex(r'$p_{turbulence}$') 
    p_pi.settex(r'$p_{\pi}$') 
    p_sig.settex(r'$p_{\sigma}$') 
    p_tau.settex(r'$p_{\tau}$')
    tracer(Run,[preal, p_turb, p_tau], phase+'pres.png')
    tracer(Run,[preal, p, p_turb, p_sig, p_tau, p_pi], phase+'pression.png')
    
    
    return res
   

def Init_parametric(FixedCase):

	FixedCase["parametric"]["uvfixe"] = [[],[]]
	FixedCase["parametric"]["ulfixe"] = [[],[]]
	FixedCase["parametric"]["urfixe"] = [[],[]]
	FixedCase["parametric"]["uvlibre"] = [[],[]]
	FixedCase["parametric"]["ullibre"] = [[],[]]
	FixedCase["parametric"]["urlibre"] = [[],[]]

	FixedCase["parametric"]["Reb"] = [[],[]]
	FixedCase["parametric"]["mu"] = [[],[]]
	FixedCase["parametric"]["cd"] = [[],[]]
	FixedCase["parametric"]["R11 WIF"] = [[],[]]
	FixedCase["parametric"]["R22 WIF"] = [[],[]]
	FixedCase["parametric"]["R33 WIF"] = [[],[]]

	FixedCase["parametric"]["R11 WIT"] = [[],[]]
	FixedCase["parametric"]["R22 WIT"] = [[],[]]
	FixedCase["parametric"]["R33 WIT"] = [[],[]]

	FixedCase["parametric"]["R11 tot"] = [[],[]]
	FixedCase["parametric"]["R22 tot"] = [[],[]]
	FixedCase["parametric"]["R33 tot"] = [[],[]]


	FixedCase["parametric"]["epsilon 11 WIF"] = [[],[]]
	FixedCase["parametric"]["redistribution 11 WIF"] = [[],[]]
	FixedCase["parametric"]["production 11 WIF"] = [[],[]]
	FixedCase["parametric"]["epsilon 22 WIF"] = [[],[]]
	FixedCase["parametric"]["redistribution 22 WIF"] = [[],[]]
	FixedCase["parametric"]["epsilon 33 WIF"] = [[],[]]
	FixedCase["parametric"]["redistribution 33 WIF"] = [[],[]]

	FixedCase["parametric"]["epsilon 11 WIT"] = [[],[]]
	FixedCase["parametric"]["redistribution 11 WIT"] = [[],[]]
	FixedCase["parametric"]["production 11 WIT"] = [[],[]]
	FixedCase["parametric"]["epsilon 22 WIT"] = [[],[]]
	FixedCase["parametric"]["redistribution 22 WIT"] = [[],[]]
	FixedCase["parametric"]["epsilon 33 WIT"] = [[],[]]
	FixedCase["parametric"]["redistribution 33 WIT"] = [[],[]]
	FixedCase["parametric"]["diffusion_3_1"] = [[],[]]
	FixedCase["parametric"]["diffusion_3_2"] = [[],[]]

	FixedCase["parametric"]["epsilon 11 tot fixe"] = [[],[]]
	FixedCase["parametric"]["redistribution 11 tot fixe"] = [[],[]]
	FixedCase["parametric"]["production 11 tot fixe"] = [[],[]]
	FixedCase["parametric"]["epsilon 22 tot fixe"] = [[],[]]
	FixedCase["parametric"]["redistribution 22 tot fixe"] = [[],[]]
	FixedCase["parametric"]["epsilon 33 tot fixe"] = [[],[]]
	FixedCase["parametric"]["redistribution 33 tot fixe"] = [[],[]]

	FixedCase["parametric"]["epsilon 11 tot libre"] = [[],[]]
	FixedCase["parametric"]["redistribution 11 tot libre"] = [[],[]]
	FixedCase["parametric"]["production 11 tot libre"] = [[],[]]
	FixedCase["parametric"]["epsilon 22 tot libre"] = [[],[]]
	FixedCase["parametric"]["redistribution 22 tot libre"] = [[],[]]
	FixedCase["parametric"]["epsilon 33 tot libre"] = [[],[]]
	FixedCase["parametric"]["redistribution 33 tot libre"] = [[],[]]
	return


def calculerReBulk(Run):
    h=1.0
    nu=Run['variable']['mu']
    rhol=Run['variable']['rho_l']
    rhov=Run['variable']['rho_v']
    Iul=Run['variable']['u_l']
    Ivuv=Run['variable']['u_v']
    I=Run['variable']['I']
    Iinv=I.inv(1e-8)
    Iv=Run['variable']['I']-1.0
    Iv=Iv*(-1.0)
    Ivinv=Iv.inv(1e-8)
    ul=Iul*Iinv
    uv=Ivuv*Ivinv
    
    debit=((ul*rhol+uv*rhov).integ())._npa[0,-1]
    reb=debit*(4*h/nu)
    vitesse=(ul.integ())._npa[0,-1]
    rebb=vitesse*(4*h/nu)
    print 'Reb=', reb
    print 'Rebb=', rebb, 'pour un ul=', vitesse
    return

def calculerEo(Run):
    h=1.0
    rhol=Run['variable']['rho_l']
    db=Run['variable']['db']
    sigma=Run['variable']['sigma']
    g=Run['variable']['g']
    print 'Eo=', rhol*db*db*g/sigma
    return
    


def PositionBulles(Run):
    pi=3.1415
    Lx=2.0*pi
    Ly=pi/2
    Lz=2.0
    Iv=((Run['variable']['I']-1.0)*(-1.0)*Lx*Ly)
    rb=Run['variable']['db']/2.0
    Iv.symetriser([1])

    #Nb=[0.2825, 0.3175, 0.3525, 0.3875, 0.4225, 0.4575, 0.4925, 0.5275, 0.5625, 0.5975, 0.6325, 0.6675, 0.7025, 0.7375, 0.7725, 0.8075, 0.8425, 0.8775, 0.9125, 0.9475, 0.9825]
    #Nb=[0.24, 0.3, 0.3475, 0.3875, 0.4225, 0.4575, 0.4925, 0.5275, 0.5625, 0.5975, 0.6325, 0.6675, 0.7025, 0.7375, 0.7725, 0.8075, 0.8425, 0.8775, 0.9125, 0.9475, 0.9825]
     #   1.76  1.7  1.6525  1.6125  1.5775  1.5425  1.5075  1.4725  1.4375  1.4025  1.3675  1.3325  1.2975  1.2625  1.2275  1.1925  1.1575  1.1225  1.0875  1.0525  1.0175
    Nb=[0.3875, 0.4225, 0.4575, 0.4925, 0.5275, 0.5625, 0.5975, 0.6325, 0.6675, 0.7025, 0.7375, 0.7725, 0.8075, 0.8425, 0.8775, 0.9125, 0.9475, 0.9825]


    Ivnew=Iv*0.0
    tracer(Run,[Iv, Ivnew], 'Ivdeb')
    for i in range(len(Nb)):

    	xb=Nb[i]
    	xbm=Lz-Nb[i]
        if xb==1.0:
		xbm=500000
        for j in range(len(Iv._ax[0,:])):
                x=Iv._ax[0,j]
        	if x>(xb-rb) and x<(xb+rb):
		        Ivnew._npa[0,j]=Ivnew._npa[0,j]+(pi*(rb*rb-(x-xb)*(x-xb)))
        	if x>(xbm-rb) and x<(xbm+rb):
		        Ivnew._npa[0,j]=Ivnew._npa[0,j]+(pi*(rb*rb-(x-xbm)*(x-xbm)))

    tracer(Run,[Iv, Ivnew], 'Iv')       
    
    return



def outDataForSrcNeptune(Run_imput):
    """ 
    print all the profils in a .txt files  
    @param a dictionnary  
    """   
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    Run=Run_imput
    for cle in Run.keys():
        for cle2 in Run[cle].keys():
            if isinstance(Run[cle][cle2], ScalarField):
                Run[cle][cle2].symetriser()
            elif isinstance(Run[cle][cle2], VectorField):
                Run[cle][cle2].symetriser()    
            elif isinstance(Run[cle][cle2], Tensor2Field):
                Run[cle][cle2].symetriser(sym_tens)                   
                
                
    L=len(Run['variable']['I']._ax[0,:])/2
    fichier = open(Run['variable']['name']+"_out_ligne.txt", 'w')
    fichier.write('# Coordonnee_z['+str(L)+'] = {')
    for j in range(L):
        if j==L-1:
            fichier.write(str(Run['variable']['I']._ax[0,j])+'};\n')
        else:
            fichier.write(str(Run['variable']['I']._ax[0,j])+', ') 
    
    for cle in Run.keys():
        for cle2 in Run[cle].keys():
            if isinstance(Run[cle][cle2], ScalarField):
                fichier.write('# '+cle+'_'+cle2+'['+str(L)+'] = {')
                for j in range(L):
                    if j==L-1:
                        fichier.write(str(Run[cle][cle2]._npa[0,j])+'};\n')
                    else:
                        fichier.write(str(Run[cle][cle2]._npa[0,j])+', ') 
            if isinstance(Run[cle][cle2], VectorField):
                for k in range(3):
                    fichier.write('# '+cle+'_'+cle2+'_'+str(k)+'['+str(L)+'] = {')
                    for j in range(L):
                        if j==L-1:
                            fichier.write(str(Run[cle][cle2]._npa[k,j])+'};\n')
                        else:
                            fichier.write(str(Run[cle][cle2]._npa[k,j])+', ') 
                    
                  
                    
            if isinstance(Run[cle][cle2], Tensor2Field):
                for k in range(3):
                    for l in range(3):
                        fichier.write('# '+cle+'_'+cle2+'_'+str(k)+str(l)+'['+str(L)+'] = {')
                        for j in range(L):
                            if j==L-1:
                                fichier.write(str(Run[cle][cle2]._npa[k,l,j])+'};\n')
                            else:
                                fichier.write(str(Run[cle][cle2]._npa[k,l,j])+', ') 
                
 
    return
    
    return
