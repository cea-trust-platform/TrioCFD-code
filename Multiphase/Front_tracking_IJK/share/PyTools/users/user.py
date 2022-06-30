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

def fluctuij(p, du, pdu):
    
    if isinstance(p, VectorField) and isinstance(du, VectorField):
        pdu_bis=p.MatProduct(du)
    if isinstance(p, ScalarField) and isinstance(du, VectorField):
        pdu_bis=p*du
    if isinstance(p, ScalarField) and isinstance(du, Tensor2Field):    
        pdu_bis=du*p
    if isinstance(p, Tensor2Field) and isinstance(du, Tensor2Field):    
        pdu_bis=p.tensorProduct(du)
        pdu_bis=pdu_bis.contract('ij', 'ikjk')
        
    res=pdu-pdu_bis
    return res

def fluctuijk(uuu, u, uu):
    Rijk=uuu
    Rij_uk=u.tensorProduct(uu)
    print Rij_uk
    uijk=u.tensorProduct(u).tensorProduct(u)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Rijk._npa[i,j,k,:]=Rijk._npa[i,j,k,:]+uijk._npa[i,j,k,:]*2.0-Rij_uk._npa[j,i,k,:]-Rij_uk._npa[i,j,k,:]-Rij_uk._npa[k,i,j,:] 
         
    #Rijk._npa[0,0,0,:]=0.0
    return Rijk    




def UsualVariable() :
    ################################################################################
    ########### PLOT DES PROFIL DE VITESSE ET INDICATRICE ##########################
    ################################################################################
    sym_tens=[1,-1,1,-1,1,-1,1,-1,1]
    sym_tens=[1,1,-1,1,1,-1,-1,-1,1]
    xlims=[0,1]
    fstat = "Stats_sphere.dt_ev"
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u=u_l+u_v
    uu_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu=uu_v+uu_l
    ufu=u.MatProduct(u)
    Rij=uu-ufu
    fstat = "Stats_deform.dt_ev"
    II=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    u2_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u2_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u2=u2_l+u2_v
    uu2_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu2_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu2=uu2_v+uu2_l
    ufu2=u2.MatProduct(u2)
    R2ij=uu2-ufu2
    fstat = "Stats_Guillaume.dt_ev"
    h=5e-3;rho=594.38;mu=6.8327e-5;Ret=180
    utau=Ret*mu/(rho*h)
    #utau=1.
    print "utau= ", utau
    I3=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I3._ax= I3._ax/h
    u3_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u3_l._ax= u3_l._ax/h
    u3_l=u3_l.product(1/utau)
    u3_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', r'$z$', r'$u_v$')
    u3_v._ax= u3_v._ax/h
    u3_v=u3_v.product(1/utau)
    u3=u3_l+u3_v
    uu3_l=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', r'$z$', r'$\langle(u\times u)\rangle)_l$')
    uu3_v=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', r'$z$', r'\langle(u\times u)\rangle)_v')
    uu3_l._ax= uu3_l._ax/h
    uu3_v._ax= uu3_v._ax/h
    uu3_l=uu3_l.product(1/(utau**2))
    uu3_v=uu3_v.product(1/(utau**2))
    uu3=uu3_v+uu3_l
    ufu3=u3.MatProduct(u3)
    R3ij=uu3-ufu3
    R3ij=R3ij.symetriser(sym_tens)
    R3ij.settex(r'WS155 (TrioCFD)')
    # Set names : 
    u_l.settex(r'$u_{liquide}$ spherical (TrioCFD)')
    u_v.settex(r'$u_{vapeur}$ spherical (TrioCFD)')
    u2_l.settex(r'$u_{liquide}$ deformable (TrioCFD)')
    u2_v.settex(r'$u_{vapeur}$ deformable (TrioCFD)')
    u3_l.settex(r'$u_{liquide}$ deformable (TrioCFD)')
    u3_l.settex(r'$u_{vapeur}$ deformable (TrioCFD)')
    I=I.add(1.0, '--')
    II=II.add(1.0, '--')
    I3=I3.add(1.0, '--')
    I_tryg=Field.LoadFromFile('I_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    I_tryg.settex(r'spherical (LT2008)')
    II_tryg=Field.LoadFromFile('I_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    II_tryg.settex(r'deformable (LT2008)')
    I.settex(r'spherical (TrioCFD)')
    II.settex(r'deformable (TrioCFD)')
    I3.settex(r'WS155 (TrioCFD)')
    I=I.symetriser()
    II=II.symetriser()
    I3=I3.symetriser()
    tracer([I, I_tryg,II, II_tryg,I3],'Void_fraction_legOff',textitle=r'Void fraction',couleur=[2,2,1], xlim=xlims, legend='off')
    tracer([I, I_tryg,II, II_tryg,I3],'Void_fraction',textitle=r'Void fraction',couleur=[2,2,1], xlim=xlims)
    
    ############################################################################################
    ############################## PLOT DES U_RMS ##############################################
    ############################################################################################
    
    # Adimensionnement par utau=47.24/2.
    Rij=Rij*(47.24/2.)*(47.24/2.)
    R2ij=R2ij*(47.24/2.)*(47.24/2.)
    #Rij=Rij.sqrt()*(47.24/2.)
    #RRij=RRij.sqrt()*(47.24/2.)
    aa=Rij
    #Rij=Rij.symetriser(sym_tens)
    R2ij=R2ij.symetriser(sym_tens)
    aa.settex(r'aa (TrioCFD)')
    Rij.settex(r'spherical (TrioCFD)')
    R2ij.settex(r'deformable (TrioCFD)')
    for i in range(2):
        for j in range(2):
            Rij.compo(i,j).settex(r'spherical (TrioCFD)')
            R2ij.compo(i,j).settex(r'deformable (TrioCFD)')
            aa.compo(i,j).settex(r'deformable (TrioCFD)')
            pass
        pass
    #
    Rij_tryg00=Field.LoadFromFile('UU_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg00=Field.LoadFromFile('UU_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Rij_tryg00=Rij_tryg00.power(2)
    R2ij_tryg00=R2ij_tryg00.power(2)
    Rij_tryg00.settex(r'spherical (LT2008)')
    R2ij_tryg00.settex(r'deformable (LT2008)')
    #tracer([Rij.compo(0,0), Rij_tryg00 , R2ij.compo(0,0), R2ij_tryg00, R3ij.compo(0,0)], 'uu',textitle=r"Axial $\overline{u'u'}$", couleur=[2,2,1], xlim=xlims)
    tracer([Rij.compo(0,0), Rij_tryg00 , R2ij.compo(0,0), R2ij_tryg00], 'uu',textitle=r"Axial $\overline{u'u'}$", couleur=[2,2], xlim=xlims, legend='off')
    #tracer([R3ij.compo(0,0)], 'uu_ws155',textitle=r"Axial $\overline{u'u'}$", couleur=[0,0,1], xlim=xlims, legend='off')
    Rij_tryg22=Field.LoadFromFile('VV_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg22=Field.LoadFromFile('VV_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    Rij_tryg22=Rij_tryg22.power(2)
    R2ij_tryg22=R2ij_tryg22.power(2)
    Rij_tryg22.settex(r'spherical (LT2008)')
    R2ij_tryg22.settex(r'deformable (LT2008)')
    tracer([Rij.compo(2,2), Rij_tryg22, R2ij.compo(2,2), R2ij_tryg22], 'vv',textitle=r"Wall-normal $\overline{v'v'}$", couleur=[2,2], xlim=xlims, legend='off')
    #tracer([Rij.compo(1,1), Rij_tryg11, R2ij.compo(1,1), R2ij_tryg11], 'ww',textitle=r"Spanwise $\overline{w'w'}$", couleur=[2,2], xlim=xlims, legend='off')
    Rij_tryg20=Field.LoadFromFile('UV_sphere_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', r'$z$', r'I')
    R2ij_tryg20=Field.LoadFromFile('UV_defo_TRYG.dt_ev',["coordonnee_K"] , ["I"],'I', 'z', r'I')
    Rij_tryg20.settex(r'spherical (LT2008)')
    R2ij_tryg20.settex(r'deformable (LT2008)')
    tracer([aa.compo(2,0)*(-1.), Rij_tryg20 , R2ij.compo(2,0)*(-1.), R2ij_tryg20], 'uv',textitle=r"Cross-correlation $\overline{u'v'}$", couleur=[2,2], xlim=xlims, legend='off')
    #tracer([Rij.compo(2,0)*(-1.), Rij_tryg20 , R2ij.compo(2,0)*(-1.), R2ij_tryg20], 'uv',textitle=r"Cross-correlation $\overline{u'v'}$", couleur=[2,2], xlim=xlims)
    #tracer([aa.compo(2,0), Rij.compo(2,0),uu.compo(2,0), ufu.compo(2,0)], 'uvbis',textitle=r"Cross-correlation $\overline{u'v'}$", xlim=xlims)
    
    return
############################################################################################
############################## CHARGEMENT DES CHAMPS #######################################
############################################################################################
def RunDiphasicCase(gravity, rho_l, rho_v, alpha_v, mu, Source, sigma, fstat):

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
    variable={'I':I,'aiN':aiN,'aii':aii, 'PaiN':PaiN, 'g':gravity, 'rho_l':rho_l, 'rho_v':rho_v, 'alpha':alpha_v,'mu': mu,'S':Source,'sigma':sigma, 'p_l':pression_l, 'p_v':pression_v}
    ############################################################################################
    ############################## PARAMETRES  #################################################
    ############################################################################################

    
    g=Field.initgravity([gravity,0,0], 'g', I)
    Iv=I.add(1.0, '--')
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


    qdm_diph=qdm_diphase('diph',u, mu, pression, Rij, S, rho, aii, rho_av, g, sigma)
    qdm_liq=qdm_diphase('liquide',u_l, mu, pression_l, Rij_l, S, rho_l, aii, rho_av, g, sigma, I)
    qdm_gaz=qdm_diphase('gaz',u_v, mu, pression_v, Rij_v, S, rho_v, aii, rho_av, g, sigma, Iv)


    ############################################################################################
    ############################## BILAN DE RIJ ################################################
    ############################################################################################
    
#     uu=u_l.tensorProduct(u_l)
#     gradu=u_l.grad()
#     Rij_l=uu_l-uu
#     Rijk_l=uuu_l-uu_l.tensorProduct(u_l)+u_l.tensorProduct(u_l).tensorProduct(u_l)*2.0
#     pu_ll=pu_l-pression_l*u_l
#     pdu_ll=pdu_l-pression_l*u_l.grad()
#     dudu0=dudxdudx_l+dudydudy_l+dudzdudz_l
#     duduijkl=gradu.tensorProduct(gradu)
#     dudu_l=dudu0-duduijkl.contract('ij', 'ikjk')
    vpainn=vPaiN-pression_l*vaiN-u_l.tensorProduct(PaiN)+pression_l*u_l.tensorProduct(aiN)
    uduain=vdvdxaiN+vdvdyaiN+vdvdzaiN
    duain=dvdxaiN+dvdyaiN+dvdzaiN
    
    ## calcul des fluctuations
    Rij_l=fluctuij(u_l, u_l, uu_l)
    pdu_ll=fluctuij(pression_l, u_l.grad(), pdu_l)
    pu_ll=fluctuij(pression_l, u_l, pu_l)
    dudu0=dudxdudx_l+dudydudy_l+dudzdudz_l
    dudu_l=fluctuij(u_l.grad(), u_l.grad(), dudu0)
    Rijk_l=fluctuijk(uuu_l, u_l, uu_l)         

    Rij_diph=Rij_diphase(u_l, Rij_l, Rijk_l, pu_ll, pdu_ll, rho_l, dudu_l, mu, vpainn, uduain, duain, vaiN, aiN, vvaiN)
    res={'qdm_diph':qdm_diph,'qdm_liq':qdm_liq,'qdm_gaz':qdm_gaz, 'Rij':Rij_diph, 'variable':variable}
    return res
############################################################################################
############################## COMPARAISON MONOPHASIQUE ####################################
############################################################################################

def RunVremanCase(fstat1, fstat2, rho, mu, Source):
    
    dvar1 = Field.getEntries(fstat1)
    dvar2 = Field.getEntries(fstat2)
    u=Field.LoadFromFile(fstat1,["y"] , ["u1", "u3", "u2"],r'$U_l$', r'$z$', r'$u_l$')
    Rij=Field.LoadFromFile(fstat2,["y"] , ['u1u1', 'u1u3', 'u1u2', 'u1u3', 'u3u3', 'u2u3', 'u1u2', 'u2u3', 'u2u2', ],'P_l', 'z', r'$P_l$')
    Rijk=Field.LoadFromFile(fstat2,["y"] , ['0','0','u1u1u2','0','0','u1u3u2','u1u1u2','u1u3u2','u1u2u2','0','0','u1u3u2','0','0','u3u3u2','u1u3u2','u3u3u2','u3u2u2','u1u1u2','u1u3u2','u1u2u2','u1u3u2','u3u3u2','u3u2u2','u1u2u2','u3u2u2','u2u2u2'],'P_l', 'z', r'$P_l$')
    pression=Field.LoadFromFile(fstat1,["y"] , ["p"],'P_l', 'z', r'$P_l$')
    pu=Field.LoadFromFile(fstat2,["y"] , ["u1p", "u3p", "u2p"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx=Field.LoadFromFile(fstat2,["y"] , ["u11u11", "u11u31", "u11u21","u11u31", "u31u31", "u21u31","u11u21", "u21u31", "u21u21"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy=Field.LoadFromFile(fstat2,["y"] , ["u13u13", "u13u33", "u13u23","u13u33", "u33u33", "u23u33","u13u23", "u23u33", "u23u23"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz=Field.LoadFromFile(fstat2,["y"] , ["u12u12", "u12u32", "u12u22","u12u32", "u32u32", "u22u32","u12u22", "u22u32", "u22u22"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu=Field.LoadFromFile(fstat2,["y"] , ["u11p","u13p","u12p","u31p","u33p","u32p","u21p","u23p","u22p"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')

    S=Field.initsource([Source,0,0], 'beta', u)
    dudu0=dudxdudx+dudydudy+dudzdudz
    u._ax[:,:]=u._ax[:,:]*mu
    Rij._ax[:,:]=Rij._ax[:,:]*mu
    Rijk._ax[:,:]=Rijk._ax[:,:]*mu
    pression._ax[:,:]=pression._ax[:,:]*mu
    pu._ax[:,:]=pu._ax[:,:]*mu
    dudxdudx._ax[:,:]=dudxdudx._ax[:,:]*mu
    dudydudy._ax[:,:]=dudydudy._ax[:,:]*mu
    dudzdudz._ax[:,:]=dudzdudz._ax[:,:]*mu
    pdu._ax[:,:]=pdu._ax[:,:]*mu
    qdm_vreman=qdm_monophase('mono', u, mu, pression, Rij, S, rho)
    Rij_vreman=Rij_monophase(u, Rij, Rijk, pu, pdu, rho, dudu0, mu)
    res={'qdm':qdm_vreman, 'Rij':Rij_vreman}
    
    return res

    
def RunMonophasicCase(fstat, rho, mu, Source):
    

    dvar = Field.getEntries(fstat)
    I=Field.LoadFromFile(fstat,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression=Field.LoadFromFile(fstat,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    uu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    aii=Field.LoadFromFile(fstat,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu=Field.LoadFromFile(fstat,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    S=Field.initsource([Source,0,0], 'beta', u)
    ## calcul des fluctuations
    Rij=fluctuij(u, u, uu)
    pdu=fluctuij(pression, u.grad(), pdu)
    pu=fluctuij(pression, u, pu)
    dudu0=dudxdudx+dudydudy+dudzdudz
    dudu0=fluctuij(u.grad(), u.grad(), dudu0)
    Rijk=fluctuijk(uuu, u, uu)
    
    qdm_mono=qdm_monophase('mono', u, mu, pression, Rij, S, rho)
    Rij_mono=Rij_monophase(u, Rij, Rijk, pu, pdu, rho, dudu0, mu)
    res={'qdm':qdm_mono, 'Rij':Rij_mono}
    return res  



#########################################################################################
###################################### MAIN #############################################
#########################################################################################

def load(fichier):
    output=open(fichier, 'rb')
    Run=pickle.load(output)
    output.close()
    clef(Run)
    return Run

def loadSansClef(fichier):
    output=open(fichier, 'rb')
    Run=pickle.load(output)
    output.close()
    return Run

def save(Var, fichier):
    output=open(fichier, 'wb')
    pickle.dump(Var,output)
    output.close()
    return

def clef(SphericalRun):
    print ''
    print 'Liste des objets'
    for cle in SphericalRun.keys():
        i=0
        try:
            for cle2 in SphericalRun[cle].keys():
                i=i+1
                if i==1:
                    print cle+'--->'+cle2
                else:
                    print '------->'+cle2 
        except:
                print '--->'+cle
           
    print ''
    return

def Mk(DeformableRun):
    
    #### chargement des donnees
    I=DeformableRun['variable']['I']
    M_l=DeformableRun['qdm_liq']['interface']
    M_v=DeformableRun['qdm_gaz']['interface']
    M_sigma_code=DeformableRun['variable']['aii']*DeformableRun['variable']['sigma']*(-1.0)
    M_sigma_code.settex('interface monofluide code')    
    p_v=DeformableRun['variable']['p_v']
    p_l=DeformableRun['variable']['p_l']
    aiN=DeformableRun['variable']['aiN']
    M_sigma_mono=DeformableRun['qdm_diph']['interface']
    
    tracer([p_v*I.grad()*(-1.0)+p_l*I.grad()], 'bla', xlim=[0,0.5*I._ax[0,-1]])
    
    DP=p_v*I.grad()*(-1.0)+p_l*I.grad()
    Mi_v=M_v+p_v*I.grad()*(-1.0)#### ATTENTION Ain change de signe en fonction de la phase
    Mi_l=M_l+p_l*I.grad()
    ###### Calculs
    
    M_sigma_bifluide=M_l+M_v
    Mi_sigma=Mi_l+Mi_v
    
    
    ###### PRINT
    DP.settex('DP')
    Mi_sigma.settex('Mi_sigma') 
    M_l.settex('M_l')
    M_v.settex('M_v')
    Mi_l.settex('Mi_l')
    Mi_v.settex('Mi_v')
    M_sigma_code.settex('M_code')
    M_sigma_bifluide.settex('M_bifluide')
    M_sigma_mono.settex('M_monofluide')
    tracer([M_l,M_v,M_sigma_code,M_sigma_bifluide, M_sigma_mono], 'M_sigma_i', [0], xlim=[0,0.5*I._ax[0,-1]])
    tracer([M_l,M_v,M_sigma_code,M_sigma_bifluide, M_sigma_mono], 'M_sigma_j', [1], xlim=[0,0.5*I._ax[0,-1]])    
    tracer([M_l,M_v,M_sigma_code,M_sigma_bifluide, M_sigma_mono], 'M_sigma_k', [2], xlim=[0,0.5*I._ax[0,-1]])
    
    tracer([Mi_l,Mi_v,Mi_sigma, M_sigma_bifluide, DP], 'Mi_sigma_i', [0], xlim=[0,0.5*I._ax[0,-1]])
    tracer([Mi_l,Mi_v,Mi_sigma, M_sigma_bifluide, DP], 'Mi_sigma_j', [1], xlim=[0,0.5*I._ax[0,-1]])    
    tracer([Mi_l,Mi_v,Mi_sigma, M_sigma_bifluide, DP], 'Mi_sigma_k', [2], xlim=[0,0.5*I._ax[0,-1]])   
    return

def tracerQdmDiff(Mono):
    tracer([Mono['qdm_diph']['interface'],Mono['qdm_diph']['viscous'],Mono['qdm_diph']['turbulence'],Mono['qdm_diph']['inertie'],Mono['qdm_diph']['pression'],Mono['qdm_diph']['source'],Mono['qdm_diph']['rhog']], 'qdm_i', [0])
    tracer([Mono['qdm_diph']['interface'],Mono['qdm_diph']['viscous'],Mono['qdm_diph']['turbulence'],Mono['qdm_diph']['inertie'],Mono['qdm_diph']['pression'],Mono['qdm_diph']['source'],Mono['qdm_diph']['rhog']], 'qdm_j', [1])
    tracer([Mono['qdm_diph']['interface'],Mono['qdm_diph']['viscous'],Mono['qdm_diph']['turbulence'],Mono['qdm_diph']['inertie'],Mono['qdm_diph']['pression'],Mono['qdm_diph']['source'],Mono['qdm_diph']['rhog']], 'qdm_k', [2])
    tracer([Mono['qdm_gaz']['interface'],Mono['qdm_gaz']['viscous'],Mono['qdm_gaz']['turbulence'],Mono['qdm_gaz']['inertie'],Mono['qdm_gaz']['pression'],Mono['qdm_gaz']['source'],Mono['qdm_gaz']['rhog']], 'qdmgaz_i', [0])
    tracer([Mono['qdm_gaz']['interface'],Mono['qdm_gaz']['viscous'],Mono['qdm_gaz']['turbulence'],Mono['qdm_gaz']['inertie'],Mono['qdm_gaz']['pression'],Mono['qdm_gaz']['source'],Mono['qdm_gaz']['rhog']], 'qdmgaz_j', [1])
    tracer([Mono['qdm_gaz']['interface'],Mono['qdm_gaz']['viscous'],Mono['qdm_gaz']['turbulence'],Mono['qdm_gaz']['inertie'],Mono['qdm_gaz']['pression'],Mono['qdm_gaz']['source'],Mono['qdm_gaz']['rhog']], 'qdmgaz_k', [2])
    tracer([Mono['qdm_liq']['interface'],Mono['qdm_liq']['viscous'],Mono['qdm_liq']['turbulence'],Mono['qdm_liq']['inertie'],Mono['qdm_liq']['pression'],Mono['qdm_liq']['source'],Mono['qdm_liq']['rhog']], 'qdmliq_i', [0])
    tracer([Mono['qdm_liq']['interface'],Mono['qdm_liq']['viscous'],Mono['qdm_liq']['turbulence'],Mono['qdm_liq']['inertie'],Mono['qdm_liq']['pression'],Mono['qdm_liq']['source'],Mono['qdm_liq']['rhog']], 'qdmliq_j', [1])
    tracer([Mono['qdm_liq']['interface'],Mono['qdm_liq']['viscous'],Mono['qdm_liq']['turbulence'],Mono['qdm_liq']['inertie'],Mono['qdm_liq']['pression'],Mono['qdm_liq']['source'],Mono['qdm_liq']['rhog']], 'qdmliq_k', [2])
        
    
    return


############################################################################################
################################## TEST ####################################################
############################################################################################


class Test(unittest.TestCase):
    
    def verif(self, a, b): 
        self.assertTrue(self.assertNumpyEqual(a._npa, b._npa))
        self.assertTrue(self.assertNumpyEqual(a._ax.all(),b._ax.all()))
        self.assertTrue(a.labelRef==b.labelRef)
        return
        
    def derivativeScalar(self):
        var=Field.initFromAnalytic('S',-5, 5, 100, 'cos', ['dsin(x)/dx'], 'x', 'blabla')
        grady = y.derivate()
        self.verif(var, grady)
        return   
    
    def derivativeVector(self):
        var=Field.initFromAnalytic('V',-5, 5, 100, 'dpoly', ['dVx/dx', 'dVy/dx', 'dVz/dx'], 'x', 'blabla')
        gradV = V.derivate()
        self.verif(var, gradV)
        return   
    
    def derivativeTensor(self):
        var=Field.initFromAnalytic('T',-5, 5, 100, 'dconst', ['dUU/dx','dUV/dx','dUW/dx','dVU/dx','dVV/dx','dVW/dx','dWU/dx','dWV/dx','dWW/dx'], 'x', 'blabla')
        gradT = T.derivate()
        self.verif(var, gradT)
        return      
     
    def addScalar(self):
        var=Field.initFromAnalytic('S',-5, 5, 100, '2sin', ['sin(x)+sin(x)'], 'x', 'blabla')
        yy = y.add(y, '+')
        self.verif(var, yy)
        return  
    
    def diffScalar(self):
        var=Field.initFromAnalytic('S',-5, 5, 100, 'zero', ['sin(x)-sin(x)'], 'x', 'blabla')
        yy = y.add(y, '-')
        self.verif(var, yy)
        return      
    
    def productScalar(self):
        var=Field.initFromAnalytic('S',-5, 5, 100, 'sin**2', ['sin(x)*sin(x)'], 'x', 'blabla')
        yy = y.product(y)    
        self.verif(var, yy)      
        return 
      
    def quotientScalar(self):
        var=Field.initFromAnalytic('S',-5, 5, 100, 'un', ['sin(x)/sin(x)'], 'x', 'blabla')
        yy =y.product(y,'/')   
        self.verif(var, yy)   
        return  
    
    def gradijk(self):
        var=Field.initFromAnalytic('T3',-5, 5, 100, 'gradsin', ['cos(x)'], 'x', 'blabla')
        var2=Field.initFromAnalytic('T4',-5, 5, 100, 'gradsin', ['cos(x)'], 'x', 'blabla')  
        y=Msin.grad()
        y2=T3.grad()
        self.verif(var, y) 
        self.verif(var2, y2)  
        return 
    
    def tensorProductijk(self):
        y0=V.tensorProduct(V)
        y1=V.tensorProduct(T)
        y2=T.tensorProduct(V)
        y3=T.tensorProduct(T)
        y4=T3.tensorProduct(V)
        var=Field.initFromAnalytic('T',-5, 5, 100, 'VV', ['cos(x)'], 'x', 'blabla')
        var1=Field.initFromAnalytic('T3',-5, 5, 100, 'VT', ['cos(x)'], 'x', 'blabla')
        var2=Field.initFromAnalytic('T3',-5, 5, 100, 'TV', ['cos(x)'], 'x', 'blabla')
        var3=Field.initFromAnalytic('T4',-5, 5, 100, 'TT', ['cos(x)'], 'x', 'blabla')
        var4=Field.initFromAnalytic('T4',-5, 5, 100, 'T3V', ['cos(x)'], 'x', 'blabla')
        self.verif(var, y0)
        self.verif(var1, y1)
        self.verif(var2, y2)
        self.verif(var3, y3)
        self.verif(var4, y4)
        return        
     
    def contractijk(self):
        y0=T4.contract('ij', 'ijkk')
        y1=T4.contract('ij', 'kijk')
        y3=T4.contract('ij', 'ikjk')
        var=Field.initFromAnalytic('T',-5, 5, 100, 'ijkk', ['cos(x)'], 'x', 'blabla')
        var1=Field.initFromAnalytic('T',-5, 5, 100, 'kijk', ['cos(x)'], 'x', 'blabla')
        var3=Field.initFromAnalytic('T',-5, 5, 100, 'ikjk', ['cos(x)'], 'x', 'blabla')
        self.verif(var, y0)
        self.verif(var1, y1)
        self.verif(var3, y3)
        return 
    
    def transpoijk(self):
        y0=T.transpo()
        var=Field.initFromAnalytic('T',-5, 5, 100, 'transpo', ['cos(x)'], 'x', 'blabla')
        self.verif(var, y0)
        return  
    
    def initfromfile(self):
        y=Field.LoadFromFile('Stats_test.dt_ev',["coordonnee_K"] , ['UUUI', 'VUUI', 'WUUI', 'UVUI', 'VVUI', 'WVUI', 'UWUI', 'VWUI', 'WWUI', 'UUVI', 'VUVI', 'WUVI', 'UVVI', 'VVVI', 'WVVI', 'UWVI', 'VWVI', 'WWVI', 'UUWI', 'VUWI', 'WUWI', 'UVWI', 'VVWI', 'WVWI', 'UWWI', 'VWWI', 'WWWI'],'UUU', 'x', r'\langle(u\times u\times u)\rangle)_v')
        var=Field.initFromAnalytic('T3',0,1,3, 'file', ['const'], 'x', 'blabla') 
#         for k in dims:
#             for j in dims:
#                 for i in dims:
#                     print y._npa[i,j,k,1]
#                     print var._npa[i,j,k,1]
#                     print ''        
        self.verif(var, y)
        return
        
#     def qdm(self):
#         vit=Field.initFromAnalytic('V',0, 2, 100, 'poiseuille', ['U', 'V', 'W'], 'x', 'blabla')
#         tauu=Field.initFromAnalytic('T',0, 2, 100, 'poiseuille', ['dU/dx','dU/dy','dU/dz','dV/dx','dV/dy','dV/dz','dW/dx','dW/dy','dW/dz'], 'x', 'blabla')
#         S=Field.initsource([-1,0,0], 'beta', vit)
#         mu=1.0      
#         tauTrans=tauu.transpo()
#         tau=tauu.add(tauTrans, '+')
#         tau=tau.product(mu, '*') 
#         divT=tau.divergence()
#         nonlin=vit.MatProduct(vit)
#         divuu=nonlin.divergence()
#         bilanqdm=divuu.add(divT, '-').add(S,'-')
#         res=bilanqdm._npa
#         res[:]=0
#         self.assertTrue(self.assertNumpyEqual(res,bilanqdm._npa, 1e-5))
         

    def assertNumpyEqual(self, a, b, delta=1e-2):
        diff = np.amax(np.abs(a-b)) 
        return diff < delta   
    
    

      

if __name__ == "__main__":
    
    
    y=Field.initFromAnalytic('S',-5, 5, 100, 'sin', ['sin(x)'], 'x', 'blabla')
    V=Field.initFromAnalytic('V',-5, 5, 100, 'poly', ['Vx', 'Vy', 'Vz'], 'x', 'blabla')
    M=Field.initFromAnalytic('V',-5, 5, 100, 'titi', ['Mx', 'My', 'Mz'], 'x', 'blabla')
    T=Field.initFromAnalytic('T',-5, 5, 100, 'const', ['UU','UV','UW','VU','VV','VW','WU','WV','WW'], 'x', 'blabla')
    T3=Field.initFromAnalytic('T3',-5, 5, 100, 'sin', ['sin'], 'x', 'blabla')
    Msin=Field.initFromAnalytic('T',-5, 5, 100, 'sin', ['sin'], 'x', 'blabla')
    T4=Field.initFromAnalytic('T4',-5, 5, 100, 'cte', ['cte'], 'x', 'blabla')
    
    suite = unittest.TestSuite()
    suite.addTest(Test('derivativeScalar'))
    suite.addTest(Test('derivativeVector'))
    suite.addTest(Test('derivativeTensor'))
    suite.addTest(Test('addScalar'))
    suite.addTest(Test('productScalar'))   
    suite.addTest(Test('diffScalar'))  
    suite.addTest(Test('quotientScalar')) 
    suite.addTest(Test('gradijk'))
    suite.addTest(Test('tensorProductijk'))
    suite.addTest(Test('contractijk'))
    suite.addTest(Test('transpoijk'))
    suite.addTest(Test('initfromfile'))
    
    unittest.TextTestRunner().run(suite)
    
    
  
    
    
    
    
    
    
    
    

