# -*- coding: utf-8 -*-
#import rlcompleter
from Tools_Rij import *
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
fstat="Stats.dt_ev"

mu = 1.39e-3 
rho = 986.51

################################################### Loading fields ################################################################################################
dvar = Field.getEntries(fstat)
I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
aiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v').contract("ijk","kji")
vaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vvaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').contract("ijk","kji")
dudxdudx_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
dudydudy_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
dudzdudz_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()

vdvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vdvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
vdvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$').transpo()
dvdxaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdyaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
dvdzaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
duaiNx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
dvaiNy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
dwaiNz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
Idu = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')

############################################# Total term ##############################################################
uddux = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IddUdxdxU", "IddVdxdxU", "IddWdxdxU", "IddUdxdxV", "IddVdxdxV", "IddWdxdxV", "IddUdxdxW", "IddVdxdxW", "IddWdxdxW"], 'Iddu', 'z', r'$IU\frac{\partial U}{\partial x}^{2}$')
udduy = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IddUdydyU", "IddVdydyU", "IddWdydyU", "IddUdydyV", "IddVdydyV", "IddWdydyV", "IddUdydyW", "IddVdydyW", "IddWdydyW"], 'Iddu', 'z', r'$IU\frac{\partial U}{\partial y}^{2}$')
udduz = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IddUdzdzU", "IddVdzdzU", "IddWdzdzU", "IddUdzdzV", "IddVdzdzV", "IddWdzdzV", "IddUdzdzW", "IddVdzdzW", "IddWdzdzW"], 'Iddu', 'z', r'$IU\frac{\partial U}{\partial z}^{2}$')
Iddudxx = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IddUdxx", "IddVdxx", "IddWdxx"], 'Iddux', 'z', r'bhooo')
Iddudyy = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IddUdyy", "IddVdyy", "IddWdyy"], 'Idduy', 'z', r'bhooo')
Iddudzz = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IddUdzz", "IddVdzz", "IddWdzz"], 'Idduz', 'z', r'bhooo')


####### Reynolds stress tensor definition #########
Rij = fluctuij(u_l,u_l,uu_l,I)
########## Useful quantities ##### 
ull=u_l*(I.inv(0.00001))
uull = uu_l*(I.inv(0.00001))
uuull = uuu_l.tensorProduct(I.inv(0.00001))
Idu_l = Idu*(I.inv(0.00001))
Idu_t = Idu.transpo()
###### Molecular diffusion (Laplacian of the reynolds stress) ####
gradRij=Rij.grad() 
ddRij0=gradRij.grad()
ddRij=ddRij0.contract('ij', 'ijkk')
##### Dissipation ######
dudu0 = dudxdudx_l + dudydudy_l + dudzdudz_l 
dull = ull.grad()
dull_t = dull.transpo()
a=dull.tensorProduct(dull)
b=a.contract("ij","ikjk")
c=Idu_t.tensorProduct(dull)
cc=c.contract("ij","kijk")
cc_t=cc.transpo() 

dudu_l = dudu0 + I*b -cc -cc_t 
dudu_l = dudu_l*(2.0)*(-1.0)
##### Interfacial terms ######
uduain=vdvdxaiN+vdvdyaiN+vdvdzaiN  
duain=dvdxaiN+dvdyaiN+dvdzaiN
aijk=ull.tensorProduct(vaiN)
ciccio=ull.tensorProduct(ull).tensorProduct(aiN)
f=aijk.grad().contract('ij', 'ijkk') ### I agree with this decomposition.. did the same ...
g=vaiN.tensorProduct(ull)
g=g.grad().contract('ij', 'ikjk')
term2= f+g

term1=vvaiN.grad()
term1=term1.contract('ij', 'ijkk')
term3=ciccio.grad()
term3=term3.contract('ij', 'ijkk')
viscinter2=(term1-term2+term3)

visc1=uduain + uduain.transpo() 
b = ull.tensorProduct(duain)
b= b + b.transpo()
c=dull.tensorProduct(vaiN)
c=c.contract('ij', 'ikjk')  
c = c + c.transpo()
visc2=b+c
e= ull.tensorProduct(ull)
de=e.grad()
deain=de.tensorProduct(aiN)
visc3=deain.contract('ij','ijkk')
viscinter1=(visc1-visc2+visc3)
viscinter = viscinter1*(-1.0) +viscinter2*(-1.0)
### Total term ####
uddux = uddux.transpo()
udduy = udduy.transpo() 
udduz = udduz.transpo() 
add = uddux + udduy + udduz
#sub_x = ull.tensorProduct(Iddudxx)
#sub_y = ull.tensorProduct(Iddudyy)
#sub_z = ull.tensorProduct(Iddudzz)
sub_x = Iddudxx
sub_y = Iddudyy
sub_z = Iddudzz
sub = sub_x + sub_y + sub_z
sub1= sub.tensorProduct(ull)
help1 = dull.grad()
help2 = help1.contract('i', 'ikk')
sub2 = help2.tensorProduct(u_l)
add2 = help2.tensorProduct(ull)*I
visc_tot = add-sub1 -sub2 + add2
visc_tot = visc_tot + visc_tot.transpo()
summ_2 = ddRij + viscinter + dudu_l
### Setting legend ###
#visc_tot.settex(r'Total $\overline{v_{i}^{,}\frac{\partial ^{2} v_{m}^{,}}{\partial x_{j}^{2}}} + \overline{v_{m}^{,}\frac{\partial ^{2} v_{i}^{,}}{\partial x_{j}^{2}}}$' )
#ddRij.settex(r'molecular diffusion  $\nu\frac{\partial^2\overline{v^,_iv^,_j}}{\partial^2 x_l}$')
#visc_tot.settex(r'$v_{i}^{,} \frac{\partial ^{2} v_{j}^{,}}{\partial x_{b}^{2}} + v_{j}^{,} \frac{\partial ^{2} v_{i}^{,}}{\partial x_{b}^{2}}$')
visc_tot.settex(r'LHS')
ddRij.settex(r'$\mathbf{D^{\mu}_{ij}}$')
#dudu_l.settex(r'dissipation $2\nu\overline{\frac{\partial v^,_i}{\partial x_l} \frac{\partial v^,_j}{\partial x_l}}$')
dudu_l.settex(r'$\epsilon_{ij}$')
#viscinter.settex(r'viscous interfacial term $\overline{\frac{\partial v^,_iv^,_jn_l\delta^I}{\partial x_l}}+\overline{\frac{\partial v^,_iv^,_j}{\partial x_l}n_l\delta^I}$')
viscinter.settex(r'$\Pi^{\mu}_{ij}$')
summ_2.settex(r'RHS')
### Plots ##
listing =[ddRij, dudu_l, viscinter, summ_2, visc_tot]
#listing=[summ_2, visc_tot]
tracer(listing, 'xx',[0])
tracer(listing, 'yy',[4])
tracer(listing, 'zz',[8])
tracer(listing, 'xz',[2])


