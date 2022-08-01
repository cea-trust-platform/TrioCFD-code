# -*- coding: utf-8 -*-

#import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sys

#%% Get Re and Pr
D=float(sys.argv[1])
Re=float(sys.argv[2])
# Pr=float(sys.argv[3])
Eo=float(sys.argv[3])

########################################################
# Calculate liquid properties
#########################################################
# Saturated water at 80 bar (liquid)
rho_liq=722.2
lambda_liq=0.55634
Cp_liq=5614.1
mu_liq=8.7766e-5

#sigma=0.015507

nu_liq=mu_liq/rho_liq
alpha_liq=lambda_liq/(rho_liq*Cp_liq)

Pr_liq=nu_liq/alpha_liq
print('Prandtl Pr: ' + str(Pr_liq))
########################################################
# Calculate vapour properties
########################################################
# Saturated water at 80 bar (vapour)
mu_vap=1.9397e-5
lambda_vap=0.067056
rho_vap=42.5
Cp_vap=1e8 #enforce zero diffusivity in the bubble

alpha_vap=lambda_vap/(rho_vap*Cp_vap)
nu_vap=mu_vap/rho_vap
print('Vapour Diffusivity: ' + str(alpha_vap))
########################################################
# Calculate the diameter to satisfy the Reynolds
########################################################
visco_ratio=mu_vap/mu_liq
print('Viscosity ratio: '+ str(visco_ratio))
if Re>5: # 5<Re<1000
    if visco_ratio<=2:
        Cd_0=48.0/Re*(1.0+2.21/np.sqrt(Re)-2.14/Re)
        Cd_2=17.0*np.power(Re,-2.0/3.0)
        Cd=(2.0-visco_ratio)/2.0*Cd_0+4.0*visco_ratio/(6.0+visco_ratio)*Cd_2
    else:
        Cd_2=17.0*np.power(Re,-2.0/3.0)
        Cd_inf=24.0/Re*(1+1.0/6.0*np.power(Re,2.0/3.0))
        Cd=4.0/(2.0+visco_ratio)*Cd_2+(visco_ratio-2.0)/(2.0+visco_ratio)*Cd_inf
else: # 1<Re<5
    Cd=8/Re*(3*visco_ratio+2)/(visco_ratio+1)*(1+0.05*(3*visco_ratio+2)/(visco_ratio+1)*Re)-0.01*(3*visco_ratio+2)/(visco_ratio+1)*Re*np.log(Re)

g=3.0/4.0*(Re*nu_liq)**2*(rho_liq/rho_vap)*Cd/(D**3)
V_D=np.pi*D**3/6.0
A_D=np.pi*D**2/4.0
v_inf=np.sqrt((2*rho_vap*V_D*g)/(rho_liq*A_D*Cd))
force_init=rho_liq*g
sigma=(rho_liq-rho_vap)*g*D**2/Eo

print('Drag coefficient Cd: ' + str(Cd))
print('Acceleration g: ' + str(g))
print('Terminal velocity v: ' + str(v_inf))
print('Surface tension: ' + str(sigma))
print('Reynolds verification: ' + str(v_inf/nu_liq*D))
########################################################
# Output the parameters into bash format 
########################################################

parameters='( "' + str(g) + '" "' + str(v_inf) + '" "' + str(Pr_liq) + '" "' + str(mu_liq) + '" "' + str(rho_liq) + '" "' + str(Cp_liq)  + '" "' + str(lambda_liq) + '" "' + str(mu_vap) + '" "' + str(rho_vap) + '" "' + str(Cp_vap)  + '" "' + str(lambda_vap) + '" "' + str(force_init) + '" "' + str(sigma) + '" )'

sys.exit(parameters)
