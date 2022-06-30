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
############################ Parameters pf the calculation ##################################
g=-9.81
rho_v=845.44
rho_l=986.51
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
vdvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","UdVdxaiNx", "VdVdxaiNx", "VdWdxaiNx","UdWdxaiNx", "VdWdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
vdvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","UdVdyaiNy", "VdVdyaiNy", "VdWdyaiNy","UdWdyaiNy", "VdWdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
vdvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","UdVdzaiNz", "VdVdzaiNz", "VdWdzaiNz","UdWdzaiNz", "VdWdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
Ivdu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IvdUdx","IvdVdx", "IvdWdx", "IvdUdy", "IvdVdy", "IvdWdy", "IvdUdz", "IvdVdz", "IvdWdz"], 'UU_l','z', r'$\nabla\cdot\u_l$')
duaiNx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
dvaiNy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
dwaiNz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
dvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
g=Field.initgravity([g,0,0], 'g', I)
########################### Quantities for the time derivate ##################################################
dvar1 = Field.getEntries(time)
u_v_1=Field.LoadFromFile(time,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_l$', r'$z$', r'$u_l$')
################################ Begin of the calculation #####################################################
uvv=u_v*(Iv.inv(0.00001))
p_v_ext = p_v_ext*(-1.0)
pll=pression_l_ext*(I.inv(0.00001))
uuv=u_v.MatProduct(u_v)
uij=uu_v-uuv
Uij=uij*rho_v
Rij=Uij
#Rij=Uij*(I.inv(0.00001))
tau_liq=duaiNx+dvaiNy+dwaiNz+ dvdxaiN+dvdyaiN+dvdzaiN ## Interfacial viscous term ###
tau_vap = tau_liq*(-1.0)

################################# Final terms in the momentum equation #####################################
tau1=Ivdu
tau2=tau1+tau1.transpo()
tau3=tau2.grad()
tau4=tau3.contract('i', 'ikk')
divT=tau4*mu

gradP=pression_v_ext.grad()
gradp=gradP*(-1.0)
divrij=Rij.grad()
divrij=divrij.contract('i','ikk')
divRij=divrij*(-1.0)
## Source term ##
av1=(rho_av-rho_v)*(-1.0) 
av2=g*av1
av3=av2*Iv
###Interfacial term###
tau_vap=tau_vap*mu
m=p_v_ext*(-1.0) + tau_vap
tau_int=tau_vap*(-1.0)
#### Inertial term ######
nonlin=uuv*rho_v
nonlin1=nonlin*(Iv.inv(0.0001))
divUU=nonlin1.grad()
divUU=divUU.contract('i','ikk')

############# Time derivative ################
nsteady=u_v_1*rho_v*(-1.0)

############ RHS ###############
ai=divT+divRij+gradp+divUU*(-1.0)+av3#+nsteady (General case: uncomment nsteady) ##### In case the source is not zero you need to add it to this term ########`
insta=ai + tau_int + p_v_ext
######### Plot #################
divUU.settex(r'advection $-\nabla \cdot \left(\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right)$') 
#divRij.settex(r'Turbulence $-\nabla \cdot \overline{(\rho\mathbf{v^,}\otimes\mathbf{v^,})}^{l}$') 
divRij.settex(r'$ \mathbf{T_{v}^{t}}$')    
#divT.settex(r'Viscous stress $\nabla \cdot [\mu \overline{(\nabla \mathbf{v}+\nabla\mathbf{v}^{t}]}^{l})$') 
divT.settex(r'$ \mathbf{T_{v}^{\mu}}$ ')  
#gradp.settex(r'Pressure $-\nabla\overline{P}^{l}$')
gradp.settex(r'$\mathbf{P_{v}}$')
#av3.settex(r'Buoyancy $(\rho-\langle\rho\rangle)\mathbf{g}$')
av3.settex(r'$\mathbf{g}$')
insta.settex(r'$\mathbf{e}$')
ai.settex(r'RHS ')
nsteady.settex(r'unsteady')
#p_l_ext.settex(r'Interfacial pressure $\langle P\nabla\chi \rangle$')
p_v_ext.settex(r'$\mathbf{M_{v}^{p}}$')
#tau_int.settex(r'Interfacial viscous stress $-\langle \tau\cdot\nabla\chi \rangle$')
tau_int.settex(r'$\mathbf{M_{v}^{\mu}}$')
#m.settex(r'Interfacial force $\langle -P\nabla\chi+\tau\cdot\nabla\chi \rangle$')
m.settex(r'$\mathbf{M_{l}}$')

 
 
stress_divRij=divRij.integ()
stress_divUU=divUU.integ()
stress_divT=divT.integ()
stress_gradp=gradp.integ()
stress_av3=av3.integ()
stress_ai=ai.integ()
stress_insta=insta.integ()
stress_m=m.integ()
stress_pext=p_v_ext.integ()     
stress_tau=tau_int.integ()

######## The terms proven to be zero in the test case (unsteady, advection, source) are notr diplayed for plots readability. Otherwise add to the following list and to the residue
tracer([divT, divRij, gradp, av3,  p_v_ext, tau_int,insta],'qdm_i', [0], markevry=[0.025]) 
tracer([divT, divRij, gradp, av3,  p_v_ext, tau_int, insta],'qdm_j', [1], markevry=[0.025])
tracer([divT, divRij, gradp, av3,  p_v_ext, tau_int, insta],'qdm_k', [2], markevry=[0.025])


tracer([m, ai], 'confront_i', [0])
tracer([m, ai], 'confront_j', [1])
tracer([m, ai], 'confront_k', [2])
 


stress_divUU.settex(r'inertia $\int^y_0\nabla .\left(\rho\overline{\mathbf{v}}\otimes\overline{\mathbf{v}}\right) dy$') 
#stress_divRij.settex(r'turbulence $\int^y_0\nabla .\left(\overline{\rho\mathbf{v^,}\otimes\mathbf{v^,}}\right) dy$')
stress_divRij.settex(r'$\int_{0}^{y}\mathbf{T}_{v}^{t}$')     
#stress_divT.settex(r'Viscous stress $\int^y_0 \nabla .(\mu \overline{(\nabla \mathbf{v}+\nabla\mathbf{v}^{t})}) dy$') 
stress_divT.settex(r'$\int_{0}^{y}\mathbf{T}_{v}^{\mu}$')  
#stress_gradp.settex(r'Pressure $\int^y_0\nabla\overline{P}dy$')
#stress_S.settex(r'source $\int^y_0\beta dy$')
#stress_av3.settex(r'Buoyancy $\int^y_0(\rho-\langle\rho\rangle)\mathbf{g}dy$')     
#stress_ai.settex(r'Interface $\int^y_0\sigma\kappa\mathbf{n}\delta^idy$')
#stress_insta.settex(r'Residual')

stress_divRij.settex(r'$\int_{0}^{y}\mathbf{T}_{v}^{t}dy$')
stress_divT.settex(r'$\int_{0}^{y}\mathbf{T}_{v}^{\mu}dy$')  
stress_gradp.settex(r'$\int_{0}^{y}\mathbf{P}_{v}dy$')
stress_av3.settex(r'$\int_{0}^{y}\mathbf{g}dy$') 
stress_insta.settex(r'$\int_{0}^{y}\mathbf{e}dy$')

#stress_ai.settex(r'RHS $-\int^y_0  \langle -P\nabla\chi+\tau\cdot\nabla\chi \rangle dy$')
#stress_m.settex(r'Interfacial contribution $-\int^y_0  \langle -P\nabla\chi+\tau\cdot\nabla\chi \rangle dy$')
#stress_pext.settex(r'Interfacial pressure contribution $\int^y_0 \langle p\nabla\chi \rangle dy$')
#stress_tau.settex(r'Interfacial viscous contibution $\int^y_0 -\langle \tau\cdot\nabla\chi \rangle dy$' )

stress_ai.settex(r'RHS' )
stress_m.settex(r'$\int_{0}^{y}M_{v}dy$')
stress_pext.settex(r'$\int_{0}^{y}\mathbf{M}_{v}^{p}dy$')
stress_tau.settex(r'$\int_{0}^{y}\mathbf{M}_{v}^{\mu}dy$' )




tracer([stress_divT, stress_divRij, stress_gradp, stress_av3, stress_pext, stress_tau,  stress_insta],'stress_qdm_i', [0], markevry=[0.025]) 
tracer([stress_divT, stress_divRij, stress_gradp, stress_av3, stress_pext, stress_tau,  stress_insta],'stress_qdm_j', [1], markevry=[0.025])
tracer([stress_divT, stress_divRij, stress_gradp, stress_av3,  stress_pext, stress_tau,  stress_insta],'stress_qdm_k', [2], markevry=[0.025])  

tracer([stress_m, stress_ai], 'stress_confront_i', [0])
tracer([stress_m, stress_ai], 'stress_confront_j', [1])
tracer([stress_m, stress_ai], 'stress_confront_k', [2])



 
