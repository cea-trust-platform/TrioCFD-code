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

clf
global compteur
compteur=0
os.system('rm -f *.png')
fstat="Stats.dt_ev"
#fstat="Stats-4.58-5.29.dt_ev"

mu = 1.39e-3 
rho = 986.51
Sinit = 9.67621403012070368e+03 
g = -9.81

############################################################## Loading Field #######################################
dvar = Field.getEntries(fstat)
I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$\alpha_l u_l$') 
uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIb", "UVI", "UWI","UVI", "VVIb", "VWI","UWI", "VWI", "WWIb"],'UU_l', 'z', r'$\alpha_l \langle(u\times u)\rangle)_l$').transpo()
aiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$n_v \delta^i$')
p_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
pression_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$\alpha_l P_l_ext$')
pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$\alpha_l P_l$')
uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUIb", "UUVIb", "UUWIb","UUVIb", "UVVIb", "UVWI","UUWIb", "UVWI", "UWWIb","UUVIb", "UVVIb", "UVWI","UVVIb", "VVVIb", "VVWIb","UVWI", "VVWIb", "VWWIb","UUWIb", "UVWI", "UWWIb","UVWI", "VVWIb", "VWWIb","UWWIb", "VWWIb", "WWWIb" ],'UUU', 'z', r'\alpha \langle(u\times u\times u)\rangle)_v').contract("ijk","kji")#### ATTT To the definition of the triple product it changes the order of interpolation/product
aii=Field.LoadFromFile(fstat,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$ n\kappa\delta^i$') 
pu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$\alpha_l up_l$')
dudxdudx_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$bad$').transpo()
dudydudy_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$bad$').transpo()
dudzdudz_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$bad$').transpo()
# p duidxj : i-ligne j-col
pdu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vPaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vvaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", \
    "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"], \
    'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').contract("ijk","kji")
p_l_ext_vaiN=Field.LoadFromFile(fstat,["coordonnee_K"], ["UP_LIQ_aiNx", "UP_LIQ_aiNy", "UP_LIQ_aiNz", "VP_LIQ_aiNx", "VP_LIQ_aiNy", "VP_LIQ_aiNz", "WP_LIQ_aiNx", "WP_LIQ_aiNy", "WP_LIQ_aiNz"], 'P_extUaiN', 'z', r'$\langle Pu\delta^i \rangle$').transpo()
pdu_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dUdy","IP_LIQ_dUdz","IP_LIQ_dVdx", "IP_LIQ_dVdy","IP_LIQ_dVdz","IP_LIQ_dWdx","IP_LIQ_dWdy","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
p_l_ext_u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UP_LIQ_I", "VP_LIQ_I", "WP_LIQ_I"],r'$Up_l$', r'$z$', r'$up_l$')
vdvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vdvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vdvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
dvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
Idu = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')


Ix=Field.LoadFromFile(fstat,["coordonnee_K"] , ["Ix"],'I', 'z', r'Ix')
xaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["xaiNx", "xaiNy", "xaiNz"],'p_l', r'$z$', r'$-x\nabla\chi_l$')
Iux=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIx", "VIx", "WIx"],r'$Up_l$', r'$z$', r'$up_l$')
vxaiN=Field.LoadFromFile(fstat,["coordonnee_K"], ["xUaiNx", "xUaiNy", "xUaiNz", "xVaiNx", "xVaiNy", "xVaiNz", "xWaiNx", "xWaiNy", "xWaiNz"], 'P_extUaiN', 'z', r'$\langle xu\delta^i \rangle$').transpo()
Ixdu = Field.LoadFromFile(fstat, ["coordonnee_K"], ["dUdxIx", "dVdxIx", "dWdxIx", "dUdyIx", "dVdyIx", "dWdyIx", "dUdzIx", "dVdzIx", "dWdzIx"], 'Ixdu', 'z', r'$\overline{\chi \nabla u}$')

########################################################################## Building of the variables ####################################################################################################

ull=u_l*(I.inv(0.00001))
dull = ull.grad()
uull = uu_l*(I.inv(0.00001))
uuull = uuu_l.tensorProduct(I.inv(0.00001))        
pression_ll_ext=pression_l_ext*(I.inv(0.00001)) 
Idu_l = Idu*(I.inv(0.00001))
Idu_t = Idu.transpo()
############### Fluctuations ###############################
# GB think It is alpha*Rij
Rij_l=fluctuij(u_l, u_l, uu_l, I)
######################### Tensors for interfacial viscous terms #################################
uduain=vdvdxaiN+vdvdyaiN+vdvdzaiN  
duain=dvdxaiN+dvdyaiN+dvdzaiN
########################### Pressure contribution ##################################
pu_int_l=p_l_ext_vaiN-pression_ll_ext*vaiN-ull.tensorProduct(p_l_ext)+pression_ll_ext*ull.tensorProduct(aiN)
pu_ll_ext=fluctuij(pression_l_ext, u_l, p_l_ext_u_l, I)
pu_ll=fluctuij(pression_l, u_l, pu_l, I)
pdu_ll_ext=fluctuij(pression_l_ext, Idu, pdu_l_ext, I)
pdu_ll = fluctuij(pression_l, Idu, pdu_l, I)

xll=Ix*(I.inv(0.00001))
xu_int_l=vxaiN-xll*vaiN-ull.tensorProduct(xaiN)+xll*ull.tensorProduct(aiN)
xu_int_l= xu_int_l*(rho*g+Sinit)
xdu_ll_ext=fluctuij(Ix, Idu, Ixdu, I)*(rho*g+Sinit)
xu = fluctuij(Ix, u_l, Iux, I)*(rho*g+Sinit)
########## Dissipation #######
dudu0=dudxdudx_l+dudydudy_l+dudzdudz_l
a=dull.tensorProduct(dull)
b=a.contract("ij","ikjk")
c=Idu_t.tensorProduct(dull)
cc=c.contract("ij","kijk")
cc_t=cc.transpo()
#dudu_l=fluctuij(Idu, Idu, dudu0, I)
dudu_l = dudu0 + I*b -cc -cc_t 
########################## Turbulent diffusion 
Rijk_l=fluctuijk(uuull, ull, uull, I).tensorProduct(I)
################################### Function for the Rij ########################################
a1=ull.tensorProduct(Rij_l) #Att! Velocita ridotta ##
b1=a1.grad()
advection=b1.contract('ij','kijk')

nonlin2=Rij_l.tensorProduct(Idu_l)
nonlin2=nonlin2.contract('ij','kijk')
nonlin3=nonlin2.transpo()
prod = (nonlin2+nonlin3) 
nonlin=advection + prod 
################################ Turbulent diffusion ###############
gradRijk=Rijk_l.grad()
nonlin4=gradRijk.contract('ij','ikjk')
###################Pressure diffusion ##########
pijk=pu_ll_ext.grad()    
p1ijk=pijk.transpo()+pijk 
xuijk= xu.grad()
xu1ijk = xuijk + xuijk.transpo()
################# Redistibution #####################
p2ijk=pdu_ll_ext.transpo()+pdu_ll_ext
xu2ijk = xdu_ll_ext + xdu_ll_ext.transpo()
################## Molecular diffusion #############
gradRij=Rij_l.grad() 
ddRij0=gradRij.grad()
ddRij=ddRij0.contract('ij', 'ijkk')
visc=ddRij*mu-dudu_l*mu*2.0 

########################## Interfacial production ############
pu_int=pu_int_l+pu_int_l.transpo()
xu_int = xu_int_l + xu_int_l.transpo()

a1=ull.tensorProduct(vaiN) #### uiuknj ####
a2=vaiN.tensorProduct(ull)
ciccio=ull.tensorProduct(ull).tensorProduct(aiN) #Probabilmente necessaria la definizione completa della velocita ull
term1=vvaiN.grad() ###### WHAT ABOUT THE THIRD ORDER TENSOR ?????? ##############
term2=a1.grad()
term2bis =a2.grad()
term3=ciccio.grad()
term1=term1.contract('ij', 'ikjk') 
term2=term2.contract('ij', 'ijkk') #### ATTT!!!! invertiti i primi due indici per come viene storato vain=ujaini
term2bis = term2bis.contract('ij','ikjk')
term3=term3.contract('ij', 'ijkk')

viscinter2=(term1 - term2 - term2bis + term3) ### Attempt to witch sign ##

visc1=uduain + uduain.transpo() 
b = ull.tensorProduct(duain) ### Changed order to have the same indices of the one implemented 
b= b + b.transpo()

c=dull.tensorProduct(vaiN)
c=c.contract('ij', 'ikjk')
c = c + c.transpo()
visc2=b+c
e=ull.tensorProduct(ull) #### The derivative if the product is considered in its whole and it is not necessary any transpo  ######
egrad= e.grad()
f = egrad.tensorProduct(aiN)   
visc3=f.contract('ij', 'ijkk')
viscinter1=(visc1-visc2+visc3) 


advection=advection*(-1.0)
production=(prod)*(-1.0)
redistribution=p2ijk*(1./rho)#+xu2ijk*(1./rho)
dissipation=dudu_l*mu*2.0*(-1.0)*(1./rho) 
diffusion_turb=nonlin4*(-1.0)
transpression=p1ijk*(-1.0)*(1./rho) #+ xu1ijk*(-1.0)*(1./rho)
diffusion_mol=ddRij*mu*(1./rho)
diffusion = diffusion_turb + diffusion_mol + transpression
interfacial_pressure=pu_int*(1./rho) #+xu_int*(1./rho)###### This sign is choerent with the formulation of Morel and the equation in the paper as they are referrede to the normal of the liquid phase, ooposed to the one provided by IJK
viscinter = (viscinter1*(-1.0) + viscinter2*(-1.0))*mu*(1./rho) 

#for i in [dissipation, transpression, redistribution, diffusion_mol, interfacial_pressure, viscinter, diffusion_turb]:
#        try:
#            i=i*(1/rho)
#        except:
#            i=i*rho.inv()

residu=(advection+production+redistribution+dissipation+transpression+diffusion_mol+diffusion_turb)*(-1.0) # e il RHS con cui confrontare il termine successivo
m=interfacial_pressure+viscinter #segno presente all'interno dell'espressione 
res_tot= residu - m

advection.settex(r'advection $\frac{DR_{i,j}}{Dt}$')
production.settex(r'production $\overline{v^,_lv^,_j}\frac{\partial\overline{v_i}}{\partial x_l}+\overline{v^,_lv^,_i}\frac{\partial\overline{v_j}}{\partial x_l}$')
redistribution.settex(r'redistribution $\overline{p^,\left(\frac{\partial v^,_j}{\partial x_i}+\frac{\partial v^,_i}{\partial x_j}\right)}$')
dissipation.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_i}{\partial x_l} \frac{\partial v^,_j}{\partial x_l}}$')
diffusion_turb.settex(r'turbulent diffusion $\frac{\partial\overline{v^,_lv^,_iv^,_j}}{\partial x_l}$')
transpression.settex(r'pressure diffusion $\frac{\partial}{\partial x_l}(\overline{p^,v^,_j}\delta_{il}+\overline{p^,v^,_i}\delta_{jl})$')
diffusion_mol.settex(r'molecular diffusion  $\nu\frac{\partial^2\overline{v^,_iv^,_j}}{\partial^2 x_l}$')
interfacial_pressure.settex(r'interfacial pressure $-\frac{1}{\rho}\overline{(P^, u_j^, n_i +P^, u_i^, n_j)\delta^i}$')
residu.settex(r'Rhs $\overline{\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}}+\overline{\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I}-\overline{ p^,v^,_jn_i\delta^I+p^,v^,_in_j\delta^I}$ ')
m.settex(r'interface $\overline{\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}}+\overline{\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I}-\overline{p^,v^,_jn_i\delta^I+p^,v^,_in_j\delta^I}$ ')
viscinter.settex(r'viscous interfacial term $\overline{\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}}+\overline{\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I}$')
res_tot.settex(r'Total residue')

listing =[advection, production, redistribution, dissipation, diffusion_turb, diffusion_mol, transpression,  interfacial_pressure, viscinter, m, residu, res_tot]
tracer(listing, 'xx',[0], markevry=[0.05])
tracer(listing, 'xy',[1], markevry=[0.05])
tracer(listing, 'xz',[2], markevry=[0.05])
tracer(listing, 'yz',[3], markevry=[0.05])
tracer(listing, 'yy',[4], markevry=[0.05])
tracer(listing, 'zz',[8], markevry=[0.05])


#Check the trace of the total interfacial production term
#And the equilibrium equation for the kinetic energy 
p1=interfacial_pressure.compo(0,0)
p2=interfacial_pressure.compo(1,1)
p3=interfacial_pressure.compo(2,2)
i1=viscinter.compo(0,0)
i2=viscinter.compo(1,1)
i3=viscinter.compo(2,2)
m1=m.compo(0,0)
m2=m.compo(1,1)
m3=m.compo(2,2)
P_1=production.compo(0,0)
P_2=production.compo(1,1)
P_3=production.compo(2,2)
trac_P=(P_1+P_2+P_3)*0.5  
r1=redistribution.compo(0,0)
r2=redistribution.compo(1,1)
r3=redistribution.compo(2,2)
red=(r1+r2+r3)*0.5
e1=dissipation.compo(0,0)
e2=dissipation.compo(1,1)
e3=dissipation.compo(2,2)
e=(e1+e2+e3)*0.5
d1=diffusion.compo(0,0)
d2=diffusion.compo(1,1)
d3=diffusion.compo(2,2)
diff_1=(d1+d2+d3)*(0.5)
t1=transpression.compo(0,0)
t2=transpression.compo(1,1)
t3=transpression.compo(2,2)
diff_p=(t1+t2+t3)*0.5
d_m_1=diffusion_mol.compo(0,0)
d_m_2=diffusion_mol.compo(1,1)
d_m_3=diffusion_mol.compo(2,2)
diff_mol=(d_m_1+d_m_2+d_m_3)*0.5
diff_tot=diff_1+diff_mol+diff_p
r1=residu.compo(0,0)
r2=residu.compo(1,1)
r3=residu.compo(2,2)
rhs=(r1+r2+r3)*0.5
r_t_1=res_tot.compo(0,0)
r_t_2=res_tot.compo(1,1)
r_t_3=res_tot.compo(2,2)
res_t=(r_t_1+r_t_2+r_t_3)*0.5
trac=(m1+m2+m3)*0.5 
trac_p=(p1+p2+p3)*0.5
trac_i=(i1+i2+i3)*0.5
    ########Plot############
trac.settex(r'$I$')
trac_p.settex(r'$I_{p}$')
trac_i.settex(r'$I_{v} $')
trac_P.settex(r'$\Pi $')
red.settex(r'$\Phi $')
e.settex(r'$\epsilon$')
diff_tot.settex(r'$D$')
rhs.settex(r'$RHS$')
res_t.settex('residu')
#tracer([trac_i, e, diff_mol], 'interface')
tracer([trac_P,red,diff_tot,e,trac,rhs,res_t], 'kinetic')




