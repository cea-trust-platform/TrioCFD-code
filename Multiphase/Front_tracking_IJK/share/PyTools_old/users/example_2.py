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

############################################################## Loading Field #######################################
dvar = Field.getEntries(fstat)
I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$\alpha_l u_l$') 
uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIb", "UVI", "UWI","UVI", "VVIb", "VWI","UWI", "VWI", "WWIb"],'UU_l', 'z', r'$\alpha_l \langle(u\times u)\rangle)_l$').transpo()
aiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$n_v \delta^i$')
p_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
pression_l_ext=Field.LoadFromFile(fstat,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$\alpha_l P_l_ext$')
pression_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$\alpha_l P_l$')
uuu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUIb", "UUVIb", "UUWIb","UUVIb", "UVVIb", "UVWI","UUWIb", "UVWI", "UWWIb","UUVIb", "UVVIb", "UVWI","UVVIb", "VVVIb", "VVWIb","UVWI", "VVWIb", "VWWIb","UUWIb", "UVWI", "UWWIb","UVWI", "VVWIb", "VWWIb","UWWIb", "VWWIb", "WWWIb" ],'UUU', 'z', r'\alpha \langle(u\times u\times u)\rangle)_v').contract("ijk","kji")
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
# That's not what we want: (wrong order)
#duaiNx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
#dvaiNy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
#dwaiNz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
duaiNx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdxaiNx","dVdxaiNx","dWdxaiNx"], 'UU_l', 'z', r'$bb$')
dvaiNy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdyaiNy","dVdyaiNy","dWdyaiNy"], 'UU_l', 'z', r'$b$')
dwaiNz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["dUdzaiNz","dVdzaiNz","dWdzaiNz"], 'UU_l', 'z', r'$bb$')
# But unused : 
del duaiNx, dvaiNy, dwaiNz
Idu = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')

############################################################# New terms for molecular difussion ( never used at the moment) #############################################################################
# uduI = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IUdUdx", "IUdVdx", "IUdWdx", "IUdUdy", "IUdVdy", "IUdWdy","IUdUdz", "IUdVdz", "IUdWdz", "IVdUdx", "IVdVdx", "IVdWdx",  "IVdUdy", "IVdVdy", "IVdWdy","IVdUdz", "IVdVdz","IVdWdz", "IWdUdx", "IWdVdx", "IWdWdx", "IWdUdy", "IWdVdy", "IWdWdy",   "IWdUdz", "IWdVdz", "IWdWdz"] ,'IUdU', 'z', r'bho' )
## T means transposed : 
# uduIT = Field.LoadFromFile(fstat, ["coordonnee_K"], ["IUdUdx", "IUdUdy", "IUdUdz", "IUdVdx", "IUdVdy", "IUdVdz", "IUdWdx", "IUdWdy", "IUdWdz", "IVdUdx", "IVdUdy", "IVdUdz", "IVdVdx", "IVdVdy", "IvdVdz", "IVdWdx", "IVdWdy", "IVdWdz","IWdUdx", "IWdUdy", "IWdUdz", "IWdVdx", "IWdVdy", "IWdVdz", "IWdWdx", "IWdWdy", "IWdWdz"], 'IUdUt','z', r'bho') 

########################################################################## Building of the variables ####################################################################################################

ull=u_l*(I.inv(0.00001))
uull = uu_l*(I.inv(0.00001))
uuull = uuu_l.tensorProduct(I.inv(0.00001))        
pression_ll_ext=pression_l_ext*(I.inv(0.00001)) 
Idu_l = Idu*(I.inv(0.00001))
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

Idu_t = Idu.transpo()
pdu_ll_ext_old=fluctuij(pression_l_ext, Idu_t, pdu_l_ext, I)
pdu_ll_old = fluctuij(pression_l, Idu_t, pdu_l, I)
# GB.begin
# Use Idu instead of the transpose Idu_t:
pdu_ll_ext=fluctuij(pression_l_ext, Idu, pdu_l_ext, I)
pdu_ll = fluctuij(pression_l, Idu, pdu_l, I)
# GB.end

# REV EN COURS ICI !!!!
############################# Dissipation 
dudu0=dudxdudx_l+dudydudy_l+dudzdudz_l
dudu_l_old=fluctuij(Idu, Idu, dudu0, I) 

# GB.begin
import pdb
from pdb import *; #pdb.set_trace()
dull = ull.grad()
dull_t = dull.transpo()
a=dull.tensorProduct(dull)
b=a.contract("ij","ikjk")

# Some bugs:
#print b.TypeField() # Sur un tensor4
# bad = dull.product(dull) # infinite Loop unprotected? -> should reject by raise Exception('unauthorized') 
# bad = 3.*I # should switch to I*3. ?
Idu_t # chi*duj/dxi
Idu # chi*dui/dxj

c=Idu_t.tensorProduct(dull)
cc=c.contract("ij","kijk")
cc_t=cc.transpo()
# dd is an alternative to cc_t :
d=dull_t.tensorProduct(Idu)
dd=d.contract("ij","kijk")
err=cc_t-dd 
if ((abs(err._npa.min())>1.e-25) or  (abs(err._npa.min())>1.e-25)):
   raise Exception("Bad")

# Trace diagonal components of err on err_out.png
# tracer([err], 'errx', [0]) 
# tracer([err], 'erry', [4]) 
# tracer([err], 'errz', [8]) 
# print err
dudu_l = dudu0 + I*b*3. -cc -cc_t
del a,b, c, d, cc, dd, cc_t
# GB.end
########################## Turbulent diffusion 
# GB think It is alpha*Rijk
Rijk_l_old=fluctuijk(uuu_l, u_l, uu_l, I) 
gradRijk_old=Rijk_l_old.grad()
nonlin4_old=gradRijk_old.contract('ij','ikjk')   
# GB.begin   
# It should be with phase-averaged variables
Rijk_l=fluctuijk(uuull, ull, uull, I).tensorProduct(I) # It's not really a tensor product as it is Tensor*Scalar
gradRijk=Rijk_l.grad()
divRijk_l=gradRijk.contract('ij','ikjk') # the first k could be in j place because the 3 firsts are permutable
nonlin4=divRijk_l
tracer([nonlin4_old], 'divRijk_l_old') 
tracer([divRijk_l], 'divRijk_l') 
# GB.end   
 
################################### Function for the Rij ########################################
# alpha is twice ...
a1_old=u_l.tensorProduct(Rij_l)
a1=ull.tensorProduct(Rij_l)
b1=a1.grad()
inertia=b1.contract('ij','kijk')
print "Why inertia is not exactly 0 but varies btw %g and %g ? "%(inertia._npa.min(), inertia._npa.max()) # Maybe numerical approx of derivatives?
#tracer([inertia], 'inertia')

# for the production, we dont want Idu_l= barre(chi dui/dxj) but instead we need d(ui_l)/dxj
nonlin2_old=Rij_l.tensorProduct(Idu_l)
nonlin2_old=nonlin2_old.contract('ij','kijk')
nonlin3_old=nonlin2_old.transpo()
prod_old = (nonlin2_old+nonlin3_old) 
# So we do:
# alpha Rij_l * dui_l/dxj
nonlin2=Rij_l.tensorProduct(dull)
nonlin2=nonlin2.contract('ij','ikjk')
nonlin3=nonlin2.transpo()
prod = (nonlin2+nonlin3) 
nonlin=inertia + prod 
tracer([prod], 'prod')
print "Why is prod so small (varies btw %g and %g) ? "%(prod._npa.min(), prod._npa.max()) # there should be the laminar part WIF

###################Pressure diffusion ##########
pijk=pu_ll_ext.grad()    
p1ijk=pijk.transpo()+pijk 
################# Redistibution #####################
p2ijk=pdu_ll_ext.transpo()+pdu_ll_ext
################## Molecular diffusion #############
gradRij=Rij_l.grad() 
ddRij0=gradRij.grad()
ddRij=ddRij0.contract('ij', 'ijkk')
# It seems visc is useless...
visc=ddRij*mu-dudu_l*mu*2.0 

########################## Interfacial production ############
pu_int=pu_int_l+pu_int_l.transpo()

a=ull.tensorProduct(vaiN)
a=a*2.
ciccio=ull.tensorProduct(ull).tensorProduct(aiN) #Probabilmente necessaria la definizione completa della velocita ull
term1=vvaiN.grad()
term2=a.grad()
term3=ciccio.grad()
term1=term1.contract('ij', 'ijkk')
term2=term2.contract('ij', 'ijkk') 
term3=term3.contract('ij', 'ijkk')
viscinter2_old=(term1-term2+term3)
del term1, term2, term3

visc1=uduain + uduain.transpo() 
b = ull.tensorProduct(duain)
b= b + b.transpo()
c=Idu_l.tensorProduct(vaiN)
c=c.contract('ij', 'kjik')
c = c + c.transpo()
visc2=b+c
e=ull.tensorProduct(aiN)
f=Idu_l.tensorProduct(e)    
f_1=f.contract('ij', 'kjik')
visc3 = f_1 + f_1.transpo()
viscinter1=(visc1-visc2+visc3)

#GB.begin
# a_ijk= ui_l barre(uj nk deltai)
aijk=ull.tensorProduct(vaiN)
ajik=aijk.contract('ijk','jik')
a=aijk+ajik
ciccio=ull.tensorProduct(ull).tensorProduct(aiN) #Probabilmente necessaria la definizione completa della velocita ull
term1=vvaiN.grad()
term2=a.grad()
term3=ciccio.grad()
term1=term1.contract('ij', 'ijkk')
term2=term2.contract('ij', 'ijkk') # Problem BUG, it should be symmetric but it's not. 

f=aijk.grad().contract('ij', 'ijkk')
g=vaiN.tensorProduct(ull)
g=g.grad().contract('ij', 'ikjk')
fg = f+g
bug=(term2+term2.transpo())*0.5-fg
print "Why do I need to symmetrize term2 to have it equal to fg ? It's a bug so I use fg which is symmetric (%g and %g) ? "%(bug._npa.min(), bug._npa.max())
# We force the use of fg:
term2=fg

term3=term3.contract('ij', 'ijkk')
# I think it's a '+' before term2...
viscinter2=(term1+term2+term3)
del term1, term2, term3, ciccio, a, aijk, ajik

visc1=uduain + uduain.transpo() 
b = ull.tensorProduct(duain)
b= b + b.transpo()
# In c, it's not Idu_l but dull
# c = d(ui_l)/dxj * barre(uk*nl deltai)   # nl -> 'l' is the 4th index.
c=dull.tensorProduct(vaiN)
c=c.contract('ij', 'ikjk')
c = c + c.transpo()
visc2=b+c
#e=ull.tensorProduct(aiN)
#f=Idu_l.tensorProduct(e)    
#f_1=f.contract('ij', 'kjik')
#visc3 = f_1 + f_1.transpo() 
e= ull.tensorProduct(ull)
de=e.grad()
deain=de.tensorProduct(aiN)
visc3=deain.contract('ij','ijkk')
# I think it's a '+' before visc2...
viscinter1=(visc1+visc2-visc3)
#GB.end

# Addition of X fields for pressure :
g=-9.81
Sinit=9.67621403012070368e+03
vb=5.235833333333333e-10
vtot=0.02*0.005*0.005
rhol=rho
rhov=845.44
alv=vb/vtot
rhom=(1-alv)*rhol+alv*rhov
deltarho=(rhol-rhov)
# The field barre(rho):
rhof = I*deltarho-rhov
fac=rhof*g+Sinit
Ifac = I*(rhol*g+Sinit)


Ix=Field.LoadFromFile(fstat,["coordonnee_K"] , ["Ix"],'I', 'z', r'Ix')
xaiN=Field.LoadFromFile(fstat,["coordonnee_K"] , ["xaiNx", "xaiNy", "xaiNz"],'p_l', r'$z$', r'$-x\nabla\chi_l$')
Iux=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIx", "VIx", "WIx"],r'$Up_l$', r'$z$', r'$up_l$')

vxaiN=Field.LoadFromFile(fstat,["coordonnee_K"], ["xUaiNx", "xUaiNy", "xUaiNz", "xVaiNx", "xVaiNy", "xVaiNz", "xWaiNx", "xWaiNy", "xWaiNz"], 'P_extUaiN', 'z', r'$\langle xu\delta^i \rangle$').transpo()

Ixfac = Ix*(rhol*g+Sinit)
Iuxfac = Iux*(rhol*g+Sinit)

xu = fluctuij(Ix, u_l, Iux, I)*(rhol*g+Sinit)
#pu_ll_ext=fluctuij(pression_l_ext, u_l, p_l_ext_u_l, I)
#pu_ll=fluctuij(pression_l, u_l, pu_l, I)
pu_ll.settex(r'pu_ll')
pu_ll_ext.settex(r'pu_ll_ext')
xu.settex(r'pux')
tracer([xu, pu_ll, pu_ll_ext], 'pux_x', [0]) 
tracer([xu, pu_ll, pu_ll_ext], 'pux_y', [1]) 
tracer([xu, pu_ll, pu_ll_ext], 'pux_z', [2]) 
print "It is significant and should go to pressure diffusion... but who cares about pressure diffusion?"

################# Redistibution #####################
# p2ijk=pdu_ll_ext.transpo()+pdu_ll_ext
pdu_ll_ext=fluctuij(pression_l_ext, Idu, pdu_l_ext, I)
pdu_ll = fluctuij(pression_l, Idu, pdu_l, I)
Ixdu = Field.LoadFromFile(fstat, ["coordonnee_K"], ["dUdxIx", "dVdxIx", "dWdxIx", "dUdyIx", "dVdyIx", "dWdyIx", "dUdzIx", "dVdzIx", "dWdzIx"], 'Ixdu', 'z', r'$\overline{\chi \nabla u}$')
xdu_ll_ext=fluctuij(Ix, Idu, Ixdu, I)*(rhol*g+Sinit)
x2ijk=xdu_ll_ext.transpo()+xdu_ll_ext
x2ijk.settex(r'x2')
pdu_ll.settex(r'pl')
pdu_ll_ext.settex(r'pl_ext')
tracer([x2ijk, pdu_ll, pdu_ll_ext], 'pdux_x', [0]) 
tracer([x2ijk, pdu_ll, pdu_ll_ext], 'pdux_y', [1]) 
tracer([x2ijk, pdu_ll, pdu_ll_ext], 'pdux_z', [2]) 


########################## Interfacial production ############
# pu_int=pu_int_l+pu_int_l.transpo()
xll=Ix*(I.inv(0.00001))
xu_int_l=vxaiN-xll*vaiN-ull.tensorProduct(xaiN)+xll*ull.tensorProduct(aiN)
xu_int_l*=(rhol*g+Sinit)

xu_int=xu_int_l+xu_int_l.transpo()
xu_int_l.settex(r'x int')
pu_int_l.settex(r'pu_int_ext')
tracer([xu_int_l, pu_int_l], 'pix_x', [0]) 
tracer([xu_int_l, pu_int_l], 'pix_y', [1]) 
tracer([xu_int_l, pu_int_l], 'pix_z', [2]) 

#totot


advection=inertia*(-1.0) # OK, it is the advection but put in RHS
production=(prod)*(-1.0) # OK, it is the Pij of Eq.5.1
redistribution=p2ijk     # OK, it is rhol*phi_ij
dissipation=dudu_l*mu*2.0*(-1.0) # OK, it is rhol*epsilon_ij
dissipation_old=dudu_l_old*mu*2.0*(-1.0) 
dissipation_old.settex(r'dissip GB $2\nu\overline{\frac{\partial v^,_i}{\partial x_l} \frac{\partial v^,_j}{\partial x_l}}$')
diffusion_turb=nonlin4*(-1.0)
transpression=p1ijk*(-1.0) 
diffusion_mol=ddRij*mu
diffusion = diffusion_turb + diffusion_mol*(1./rho) + transpression*(1./rho) # OK, it is Dij of Eq.5.1
#interfacial_pressure=pu_int
#viscinter = (viscinter1*(-1.0) +viscinter2*(-1.0))*mu
# I hate signs!!!
interfacial_pressure=pu_int*(-1.0)
viscinter = (viscinter1*(-1.0) +viscinter2*(+1.0))*mu
production_interf = (interfacial_pressure+viscinter)*(1./rho) # OK, it is Pi_ij of Eq.5.1

###################################
# this loop only changes "i", not the fields... 
# it acts on copies... 
###################################
# Diffusion_turb shouldnot be in the list:
# for i in [dissipation, dissipation_GB, transpression, redistribution, diffusion_mol, interfacial_pressure, viscinter, diffusion_turb]:
if 0: 
   for i in [dissipation, dissipation_old, transpression, redistribution, diffusion_mol, interfacial_pressure, viscinter]:
        try:
            i=i*(1./rho)
        except:
            raise Exception("rho is not a field, so you can't call rho.inv()")
            i=i*rho.inv()
else:
   dissipation *=(1./rho)
   dissipation_old *=(1./rho)
   transpression *=(1./rho)
   redistribution *=(1./rho)
   diffusion_mol *=(1./rho)
   interfacial_pressure *=(1./rho)
   viscinter *=(1./rho)
   pass

# residu=(advection+production+redistribution+dissipation+transpression+diffusion_mol)*(-1.0) # e il RHS con cui confrontare il termine successivo
# Why was diffusion_turb missing?
residu=(advection+production+redistribution+dissipation+diffusion_turb+transpression+diffusion_mol)*(-1.0) # e il RHS con cui confrontare il termine successivo
m=interfacial_pressure+viscinter #segno presente all'interno dell'espressione # OK, it is Pi_ij of Eq.5.1
mmax=max(abs((production_interf - m )._npa.min()), abs((production_interf - m )._npa.max()))
if (mmax > 1.e-15):
   raise Exception("production_interf and m should be the same and the difference is of %g"%mmax)

#I Think it's a "+" because we should have in a perfect world: residu=-production_interf
res_tot= residu - m # OK, residu and res_tot are defined on LHS of Eq.5.1. res_tot should be 0. And residue is the old way of evaluating the interfacial production.

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
tracer(listing, 'yy',[4], markevry=[0.05])
tracer(listing, 'zz',[8], markevry=[0.05])

# Here I am... 
toto
exit(0)

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
tracer([trac_i, e, diff_mol], 'interface')
#tracer([trac_P,red,diff_tot,e,trac,rhs,res_t], 'kinetic')






