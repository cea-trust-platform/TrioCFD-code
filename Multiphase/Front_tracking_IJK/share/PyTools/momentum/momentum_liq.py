# -*- coding: utf-8 -*-
import os, sys, readline, rlcompleter
readline.parse_and_bind('tab:complete')
import DNSTools as dtool
from Tools import *
import unittest
import pickle
from matplotlib.pyplot import clf
#from scipy.sparse.linalg.isolve.minres import Ainfo


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
############################ Parameters of the calculation (Current settings for the single rising bubble) ##################################
fstat = "Stats.dt_ev" 
time = "Unsteady.dt_ev"
simu = dtool.Simu(jdd='DNS.data', repr_file='diph.sauv', Ret=0.)
alpha_l = (1.0-simu.alv)
    
################# Loading Averaged fields #####################################################
dvar = Field.getEntries(fstat)
I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
Iv=(I-1.0)*(-1.0)
pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
pression_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
pression_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
pression_v_ext=Field.LoadFromFile(fstat,["coordonnee_K"], ["P_VAP_Iv"], 'P_VAP_Iv', 'z', r'$P_v_ext$')
p_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
p_v_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],'p_v', r'$z$', r'$p_v\nabla\chi_v$')
u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
aii=Field.LoadFromFile(fstat,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
aiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
vdvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","UdVdxaiNx", "VdVdxaiNx", "VdWdxaiNx","UdWdxaiNx", "VdWdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vdvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","UdVdyaiNy", "VdVdyaiNy", "VdWdyaiNy","UdWdyaiNy", "VdWdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vdvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","UdVdzaiNz", "VdVdzaiNz", "VdWdzaiNz","UdWdzaiNz", "VdWdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
IdU =  Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'UU_l', 'z', r'$ \nabla\cdot u_l$')
duaiNx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
dvaiNy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
dwaiNz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
dvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
g=Field.initgravity([simu.g,0,0], 'g', I)
########################### Field for the time derivate ##################################################
if os.path.isfile(time):
   dvar1 = Field.getEntries(time)
   u_l_1=Field.LoadFromFile(time,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
   ############# Time derivative ################
   nsteady=u_l_1*simu.rhol*(-1.0) 
else:
   print "Unsteady term is assumed nul as file is missing!"
   nsteady=u_l*(0.0)
   pass

################################ Begin of the calculation #####################################################
ull=u_l*(I.inv(0.00001))
pll=pression_l_ext*(I.inv(0.00001))
uull=uu_l*(I.inv(0.00001))
uul=ull.MatProduct(ull)
uij=uull-uul
Uij=uij*simu.rhol
Rij=Uij
#Rij=Uij*(I.inv(0.00001))
tau_liq=duaiNx+dvdxaiN+dvaiNy+dvdyaiN+dwaiNz+dvdzaiN
################################# Final terms in the momentum equation #####################################
tau = IdU + IdU.transpo()
tau3=tau.grad()
tau4=tau3.contract('i', 'ikk')
divT=tau4*simu.mul
tau3.settex(r'tau3' )
tracer([tau3], 'tau3o', [14,20], markevry=[0.0002])
tracer([tau4], 'tau4', [0], markevry=[0.0002])
tracer([divT], 'divT', [0], markevry=[0.0002])
tracer([ull], 'ullx', [0], markevry=[0.0002])
tracer([ull], 'ully', [1], markevry=[0.0002])
tracer([ull], 'ullz', [2], markevry=[0.0002])
# toto

gradP=pression_l_ext.grad()
gradp=gradP*(-1.0)
divrij=Rij.grad()
divrij=divrij.contract('i','ikk')
divRij=divrij*(-1.0)
## Source term ##
St_value=0.
rhol_g=simu.g*simu.rhol
S_theo = -simu.g*simu.rhom
if os.path.isfile("mean_acc.txt"):
   mean_acc = np.loadtxt("mean_acc.txt")[2]
   St_value=mean_acc-S_theo
   print "rho*g+S should be %g but has a true mean=%g.\n Consequently, the difference is St=%g"%(S_theo, mean_acc, St_value)
   pass
else: 
   raise Exception("mean_acc.txt missing, you should generate it from *_acceleration.out")

beta_field=Field.initgravity([rhol_g+S_theo,0,0], 'rhol_g+S_theo', I)
av3=I*beta_field

St=Field.initgravity([St_value,0,0], 'St', I)
#av3_b = I*((-1.0)*(g)*simu.alv*(simu.rhol-simu.rhov))+I*(St)
av3b = I*((g)*(simu.rhom-simu.rhol)+(St))
av3b = av3b * (-1.0)
av3b.settex(r'$\mathbf{S_g}$ (reel)')
av3n = I*(St)
av3n.settex(r'$\mathbf{\varepsilon}_n$')

# import pdb; pdb.set_trace()

###Interfacial term###
tau_liq=tau_liq*simu.mul
m=p_l_ext*(-1.0) + tau_liq
tau_int=tau_liq*(-1.0)

#### Inertial term ######
nonlin=uul*simu.rhol
nonlin1=nonlin*(I.inv(0.0001))
divUU=nonlin1.grad()
divUU=divUU.contract('i','ikk')

############ RHS ###############
ai=divT+divRij+gradp+divUU*(-1.0)+av3+av3n+nsteady
insta=ai + tau_int + p_l_ext
######### Plot #################
divUU.settex(r'advection $-\nabla \cdot \left(\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right)$') 
#divRij.settex(r'Turbulence $-\nabla \cdot \overline{(\rho\mathbf{v^,}\otimes\mathbf{v^,})}^{l}$') 
divRij.settex(r'$ \mathbf{T_{l}^{t}}$')    
#divT.settex(r'Viscous stress $\nabla \cdot [\mu \overline{(\nabla \mathbf{v}+\nabla\mathbf{v}^{t}]}^{l})$') 
divT.settex(r'$ \mathbf{T_{l}^{\mu}}$ ')  
#gradp.settex(r'Pressure $-\nabla\overline{P}^{l}$')
gradp.settex(r'$\mathbf{P_{l}}$')
#av3.settex(r'Buoyancy $(\rho-\langle\rho\rangle)\mathbf{g}$')
av3.settex(r'$\mathbf{S_g}$ (theo)')
insta.settex(r'$\mathbf{e}$')
ai.settex(r'RHS ')
nsteady.settex(r'unsteady')
#p_l_ext.settex(r'Interfacial pressure $\langle P\nabla\chi \rangle$')
p_l_ext.settex(r'$\mathbf{M_{l}^{p}}$')
#tau_int.settex(r'Interfacial viscous stress $-\langle \tau\cdot\nabla\chi \rangle$')
tau_int.settex(r'$\mathbf{M_{l}^{\mu}}$')
#m.settex(r'Interfacial force $\langle -P\nabla\chi+\tau\cdot\nabla\chi \rangle$')
m.settex(r'$\mathbf{M_{l}}$')

tracer([divT, divRij, gradp, av3, av3b, av3n, p_l_ext, tau_int, insta, nsteady],'qdm_i', [0],markevry=[0.025]) 
tracer([divT, divRij, gradp, av3, av3b,  p_l_ext, tau_int, insta, nsteady],'qdm_j', [1],markevry=[0.025])
tracer([divT, divRij, gradp, av3, av3b,  p_l_ext, tau_int, insta, nsteady],'qdm_k', [2],markevry=[0.025])

tracer([m, ai], 'confront_i', [0], markevry=[0.0002])
tracer([m, ai], 'confront_j', [1])
tracer([m, ai], 'confront_k', [2])

tracer([av3, av3b], 'compa_i', [0], markevry=[0.0002])

stress_divRij=divRij.integ()
stress_divUU=divUU.integ()
stress_divT=divT.integ()
stress_gradp=gradp.integ()
stress_av3=av3.integ()
stress_ai=ai.integ()
stress_insta=insta.integ()
stress_m=m.integ()
stress_pext=p_l_ext.integ()     
stress_tau=tau_int.integ()
stress_nsteady=nsteady.integ()
stress_divUU.settex(r'inertia $\int^y_0\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right) dy$') 
#stress_divRij.settex(r'turbulence $\int^y_0\nabla .\left(\overline{\rho\mathbf{v^,}\otimes\mathbf{v^,}}\right) dy$')
stress_divRij.settex(r'$\int_{0}^{y}\mathbf{T}_{l}^{t} dy $')     
#stress_divT.settex(r'Viscous stress $\int^y_0 \nabla .(\mu \overline{(\nabla \mathbf{v}+\nabla\mathbf{v}^{t})}) dy$') 
stress_divT.settex(r'$\int_{0}^{y}\mathbf{T}_{l}^{\mu} dy$')  
#stress_gradp.settex(r'Pressure $\int^y_0\nabla\overline{P}dy$')
stress_gradp.settex(r'$\int_{0}^{y}\mathbf{P}_{l} dy$')
#stress_S.settex(r'source $\int^y_0\beta dy$')
stress_av3.settex(r'$\int_{0}^{y}\mathbf{g} dy$')     
#stress_ai.settex(r'Interface $\int^y_0\sigma\kappa\mathbf{n}\delta^idy$')
stress_insta.settex(r'$\int_{0}^{y}\mathbf{e} dy$')
stress_nsteady.settex(r'$\int_{0}^{y}\partial_t \rhol \overline{\mathbf{u}}^l dy$')

#stress_ai.settex(r'RHS $-\int^y_0  \langle -P\nabla\chi+\tau\cdot\nabla\chi \rangle dy$')
#stress_m.settex(r'Interfacial contribution $-\int^y_0  \langle -P\nabla\chi+\tau\cdot\nabla\chi \rangle dy$')
#stress_pext.settex(r'Interfacial pressure contribution $\int^y_0 \langle p\nabla\chi \rangle dy$')
#stress_tau.settex(r'Interfacial viscous contibution $\int^y_0 -\langle \tau\cdot\nabla\chi \rangle dy$' )

stress_ai.settex(r'RHS' )
stress_m.settex(r'$\int_{0}^{y} M_{l} dy$')
stress_pext.settex(r'$\int_{0}^{y} \mathbf{M}_{l}^{p} dy$')
stress_tau.settex(r'$\int_{0}^{y} \mathbf{M}_{l}^{\mu} dy$' )

tracer([stress_divT, stress_divRij, stress_gradp, stress_av3, stress_pext, stress_tau, stress_insta, stress_nsteady],'stress_qdm_i', [0], markevry=[0.025]) 
tracer([stress_divT, stress_divRij, stress_gradp, stress_av3, stress_pext, stress_tau, stress_insta, stress_nsteady],'stress_qdm_j', [1], markevry=[0.025])
tracer([stress_divT, stress_divRij, stress_gradp, stress_av3,  stress_pext, stress_tau, stress_insta, stress_nsteady],'stress_qdm_k', [2], markevry=[0.025])  

tracer([stress_m, stress_ai], 'stress_confront_i', [0], markevry=[0.04])
tracer([stress_m, stress_ai], 'stress_confront_j', [1])
tracer([stress_m, stress_ai], 'stress_confront_k', [2])

 
