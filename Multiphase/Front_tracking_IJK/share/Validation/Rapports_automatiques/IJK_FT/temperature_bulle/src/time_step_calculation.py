# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:11:47 2020

@author: mg263705
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

nbelem = float(sys.argv[1]) 
ldomain = float(sys.argv[2]) 
fo = float(sys.argv[3])
 
rho_0 = float(sys.argv[4]) 
Cp_0 = float(sys.argv[5]) 
lambda_0 = float(sys.argv[6]) 
 
rho_1 = float(sys.argv[7])
Cp_1 = float(sys.argv[8]) 
lambda_1 = float(sys.argv[9])

tf_scope=float(sys.argv[10])

#Diffusion coefficient
alpha_0 = lambda_0/(rho_0*Cp_0)
alpha_1 = lambda_1/(rho_1*Cp_1)

dt=1./6.*fo*(np.power(ldomain/nbelem,2))/(max(alpha_0,alpha_1))
    
nb_dt=np.floor(tf_scope/dt)
nb_dt=str(int(nb_dt + (100-np.mod(nb_dt,100)))) 
centroid=str(ldomain/2)

tab_exit= '( "' + nb_dt + '" "' + centroid + '" )'

sys.exit(tab_exit)
