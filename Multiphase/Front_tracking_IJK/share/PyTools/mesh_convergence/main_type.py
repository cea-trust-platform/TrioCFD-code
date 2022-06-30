# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, "/home/ab256925/PUBLIC/Pr_Alan/")

from Tools import *
from user import *
import readline
import unittest
import os
import pickle
from matplotlib.pyplot import clf
readline.parse_and_bind('tab:complete')
clf
global compteur
compteur=0
os.system('rm -f *.png')
print ''


##################################################################################################################################
########################################################## CALCULS ############################################################
##################################################################################################################################
if 1:
    #gravity=-9.81  
    #rho_l=986.51
    #rho_v=845.44
    #alpha_v=0.0010472
    #mul=1.39e-3
    #muv=1.39e-3
    #Source=0.
    #sigma=0.00163  #sigma=0.00512    
    #fstat1 = "Stats_coarse.dt_ev"   
    #time1 = "Nsteady_coarse.dt_ev"
    #fstat1 = "Stats_inter.dt_ev"   
    #time1 = "Nsteady_inter.dt_ev"
    DeformableRun=Residue(fstat1,fstat2,time1,time2,gravity, rho_l, rho_v, alpha_v, mul, muv, Source, sigma)
    #save(DeformableRun, 'data_DeformableRun')
     


##################################################################################################################################
########################################################## TRAITEMENT ############################################################
##################################################################################################################################


DeformableRun=load('data_DeformableRun')
DeformableRun['variable']['fstat'] = 'Stats_deform.dt_ev'

if 0:
	
	DeformableRun['variable']['name']='Retau180Eo3.74'
	#outDataForStat(DeformableRun)
        DeformableRun['variable']['I']=(DeformableRun['variable']['I']-1.)*(-1.)
	outDataForAutovnv(DeformableRun)

### Test a priori ###
if 1:
	DeformableRun['variable']['db'] = 0.306
	#calculerEo(DeformableRun)
	#calculerReBulk(DeformableRun)
	DeformableRun['variable']['markevery'] = [0.03]
	DeformableRun['variable']['name'] = 'Deform180'
	DeformableRun['variable']['fstat'] = "Stats.dt_ev"
	fstat = 'profil_QdM.dat'
        Check(DeformableRun)
	#CalculerKappaEtoile(DeformableRun)
        #CalculerKappaNprime(DeformableRun)
	#MlVsForcesApriori(DeformableRun)
 	#EquationVrDNS_old(DeformableRun)
	#EquationVrNeptune(fstat, DeformableRun)

