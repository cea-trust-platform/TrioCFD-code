# -*- coding: utf-8 -*-
#import rlcompleter
from Tools import *
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

######### Settings #################### 
g=-9.81 
rho_l=986.51
rho_v=845.44 
alpha_v=0.0010472
mu=1.39e-3
sigma=0.00163 
fstat = "Stats.dt_ev" 
alpha_l = (1.0-alpha_v)
rho_av = alpha_v*rho_v + alpha_l*rho_l
time = "Nsteady.dt_ev"    
################# Loading Averaged fields #####################################################
dvar = Field.getEntries(fstat)
I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
RHS = Field.LoadFromFile(fstat,["coordonnee_K"], ["IdUdx","IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
err = Field.LoadFromFile(fstat, ["coordonnee_K"], ["UaiNx", "VaiNx", "WaiNx", "UaiNy", "VaiNy", "WaiNy", "UaiNz", "VaiNz", "WaiNz"],'UU_l', 'z', r'$\nabla\cdot\u_l$' )
######### Actual Tau #######################
tau1=u_l.grad()
err = err*(-1.0)
residue= tau1  +err + RHS*(-1.0)
RHS.settex(r'$ \overline{\chi \nabla u}$')
tau1.settex(r'$ \nabla \overline{ \chi u}$')
err.settex(r'$ \overline{u \nabla \chi}$')
residue.settex(r'residue')

tracer([RHS,tau1,err, residue],'check_i', [0]) 
tracer([RHS,tau1,err, residue],'check_j', [1]) 
tracer([RHS,tau1,err, residue],'check_k', [2]) 


