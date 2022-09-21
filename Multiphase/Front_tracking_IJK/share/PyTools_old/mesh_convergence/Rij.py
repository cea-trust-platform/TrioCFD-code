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
fstatc="Stats_coarse.dt_ev"
fstatf="Stats_fin.dt_ev"
fstati="Stats_int_3.dt_ev"
fstat1="Stats_int_1.dt_ev"
fstat2="Stats_int_2.dt_ev"
mu = 1.39e-3 
rho = 986.51
def Residue2(fstac,fstati,fstaf, fstat1, fstat2,rho,mu):
########################################################## Loading Fields #############################################################################################################################
################################################################# First mesh, indicated by c ##########################################################################################################
    dvar_c = Field.getEntries(fstatc)
    I_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    u_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    uu_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    aiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    p_l_ext_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    pression_l_ext_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$P_l_ext$')
    uu_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    aii_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudxdudx_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdUdy", "IdUdxdUdx","IdUdxdUdy", "IdUdydUdy", "IdUdydUdz","IdUdxdUdz", "IdUdydUdz", "IdUdzdUdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudydudy_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdVdxdVdx", "IdVdxdVdy", "IdVdxdVdz","IdVdxdVdy", "IdVdydVdy", "IdVdydVdz","IdVdxdVdz", "IdVdydVdz", "IdVdzdVdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudzdudz_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdWdxdWdx", "IdWdxdWdy", "IdWdxdWdz","IdWdxdWdy", "IdWdydWdy", "IdWdydWdz","IdWdxdWdz", "IdWdydWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    p_l_ext_vaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"], ["UP_LIQ_aiNx", "UP_LIQ_aiNy", "UP_LIQ_aiNz", "VP_LIQ_aiNx", "VP_LIQ_aiNy", "VP_LIQ_aiNz", "WP_LIQ_aiNx", "WP_LIQ_aiNy", "WP_LIQ_aiNz"], 'P_extUaiN', 'z', r'$\langle Pu\delta^i \rangle$')
#pdu_l_ext=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dVdx","IP_LIQ_dWdx","IP_LIQ_dUdy", "IP_LIQ_dVdy","IP_LIQ_dWdy","IP_LIQ_dUdz","IP_LIQ_dVdz","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    pdu_l_ext_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dUdy","IP_LIQ_dUdz","IP_LIQ_dVdx", "IP_LIQ_dVdy","IP_LIQ_dVdz","IP_LIQ_dWdx","IP_LIQ_dWdy","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    p_l_ext_u_l_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UP_LIQ_I", "VP_LIQ_I", "WP_LIQ_I"],r'$Up_l$', r'$z$', r'$up_l$')

    vdvdxaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    duaiNx_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_c=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    Idu_c = Field.LoadFromFile(fstatc, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
####################################################################### Constructuion of the terms #################################################################################################
    ull_c=u_l_c*(I_c.inv(0.00001))
    uull_c = uu_l_c*(I_c.inv(0.00001))
 
    Rij_l_c=fluctuij(u_l_c, u_l_c, uu_l_c, I_c)
    uduain_c=vdvdxaiN_c+vdvdyaiN_c+vdvdzaiN_c  #Tensore che ha a che fare con il secondo termine viscoso all'interfaccia
    duain_c=dvdxaiN_c+dvdyaiN_c+dvdzaiN_c
    pression_ll_ext_c=pression_l_ext_c*(I_c.inv(0.00001)) 
    pu_int_l_c=p_l_ext_vaiN_c-pression_ll_ext_c*vaiN_c-ull_c.tensorProduct(p_l_ext_c)+pression_ll_ext_c*ull_c.tensorProduct(aiN_c)
    Idu_t_c = Idu_c.transpo()
    pdu_ll_ext_c=fluctuij(pression_l_ext_c, Idu_t_c, pdu_l_ext_c, I_c)
    pdu_ll_c = fluctuij(pression_l_c, Idu_t_c, pdu_l_c, I_c)
    pu_ll_ext_c=fluctuij(pression_l_ext_c, u_l_c, p_l_ext_u_l_c, I_c)
    pu_ll_c=fluctuij(pression_l_c, u_l_c, pu_l_c, I_c)
    dudu0_c=dudxdudx_l_c+dudydudy_l_c+dudzdudz_l_c ##This is a tensor

    dudu_l_c=fluctuij(Idu_c, Idu_c, dudu0_c, I_c) #Fluttuazione del gradiente di velocita 
    Rijk_l_c=fluctuijk(uuu_l_c, u_l_c, uu_l_c) #Prodotto delle velocita all'interno del termine di diffusione 

    Rij_c=Rij(u_l_c, ull_c, Rij_l_c, Rijk_l_c, pu_ll_ext_c, pdu_ll_ext_c, rho, dudu0_c, mu, I_c, uduain_c, duain_c, vaiN_c, aiN_c, vvaiN_c, pu_int_l_c, Idu_c)

################################################################# first intermediate mesh ########################################################################################################
    dvar_1 = Field.getEntries(fstat1)
    I_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    u_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    uu_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    aiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    p_l_ext_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    pression_l_ext_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$P_l_ext$')
    uu_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    aii_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudxdudx_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdUdy", "IdUdxdUdx","IdUdxdUdy", "IdUdydUdy", "IdUdydUdz","IdUdxdUdz", "IdUdydUdz", "IdUdzdUdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudydudy_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdVdxdVdx", "IdVdxdVdy", "IdVdxdVdz","IdVdxdVdy", "IdVdydVdy", "IdVdydVdz","IdVdxdVdz", "IdVdydVdz", "IdVdzdVdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudzdudz_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdWdxdWdx", "IdWdxdWdy", "IdWdxdWdz","IdWdxdWdy", "IdWdydWdy", "IdWdydWdz","IdWdxdWdz", "IdWdydWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    p_l_ext_vaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"], ["UP_LIQ_aiNx", "UP_LIQ_aiNy", "UP_LIQ_aiNz", "VP_LIQ_aiNx", "VP_LIQ_aiNy", "VP_LIQ_aiNz", "WP_LIQ_aiNx", "WP_LIQ_aiNy", "WP_LIQ_aiNz"], 'P_extUaiN', 'z', r'$\langle Pu\delta^i \rangle$')
#pdu_l_ext=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dVdx","IP_LIQ_dWdx","IP_LIQ_dUdy", "IP_LIQ_dVdy","IP_LIQ_dWdy","IP_LIQ_dUdz","IP_LIQ_dVdz","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    pdu_l_ext_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dUdy","IP_LIQ_dUdz","IP_LIQ_dVdx", "IP_LIQ_dVdy","IP_LIQ_dVdz","IP_LIQ_dWdx","IP_LIQ_dWdy","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    p_l_ext_u_l_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UP_LIQ_I", "VP_LIQ_I", "WP_LIQ_I"],r'$Up_l$', r'$z$', r'$up_l$')

    vdvdxaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    duaiNx_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_1=Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    Idu_1 = Field.LoadFromFile(fstat1, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
####################################################################### Constructuion of the terms #################################################################################################
    ull_1=u_l_1*(I_1.inv(0.00001))
    uull_1 = uu_l_1*(I_1.inv(0.00001))
 
    Rij_l_1=fluctuij(u_l_1, u_l_1, uu_l_1, I_1)
    uduain_1=vdvdxaiN_1+vdvdyaiN_1+vdvdzaiN_1  #Tensore che ha a che fare con il secondo termine viscoso all'interfaccia
    duain_1=dvdxaiN_1+dvdyaiN_1+dvdzaiN_1
    pression_ll_ext_1=pression_l_ext_1*(I_1.inv(0.00001)) 
    pu_int_l_1=p_l_ext_vaiN_1-pression_ll_ext_1*vaiN_1-ull_1.tensorProduct(p_l_ext_1)+pression_ll_ext_1*ull_1.tensorProduct(aiN_1)
    Idu_t_1 = Idu_1.transpo()
    pdu_ll_ext_1=fluctuij(pression_l_ext_1, Idu_t_1, pdu_l_ext_1, I_1)
    pdu_ll_1 = fluctuij(pression_l_1, Idu_t_1, pdu_l_1, I_1)
    pu_ll_ext_1=fluctuij(pression_l_ext_c, u_l_1, p_l_ext_u_l_1, I_1)
    pu_ll_1=fluctuij(pression_l_1, u_l_1, pu_l_1, I_1)
    dudu0_1=dudxdudx_l_1+dudydudy_l_1+dudzdudz_l_1 ##This is a tensor

    dudu_l_1=fluctuij(Idu_1, Idu_1, dudu0_1, I_1) #Fluttuazione del gradiente di velocita 
    Rijk_l_1=fluctuijk(uuu_l_1, u_l_1, uu_l_1) #Prodotto delle velocita all'interno del termine di diffusione 

    Rij_1=Rij(u_l_1, ull_1, Rij_l_1, Rijk_l_1, pu_ll_ext_1, pdu_ll_ext_1, rho, dudu0_1, mu, I_1, uduain_1, duain_1, vaiN_1, aiN_1, vvaiN_1, pu_int_l_1, Idu_1)
################################################################# Second mesh, indicated by i #######################################################################################################

    dvar_i = Field.getEntries(fstati)
    I_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    u_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    uu_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    aiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    p_l_ext_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    pression_l_ext_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$P_l_ext$')
    uu_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    aii_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudxdudx_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdUdy", "IdUdxdUdx","IdUdxdUdy", "IdUdydUdy", "IdUdydUdz","IdUdxdUdz", "IdUdydUdz", "IdUdzdUdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudydudy_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdVdxdVdx", "IdVdxdVdy", "IdVdxdVdz","IdVdxdVdy", "IdVdydVdy", "IdVdydVdz","IdVdxdVdz", "IdVdydVdz", "IdVdzdVdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudzdudz_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdWdxdWdx", "IdWdxdWdy", "IdWdxdWdz","IdWdxdWdy", "IdWdydWdy", "IdWdydWdz","IdWdxdWdz", "IdWdydWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    p_l_ext_vaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"], ["UP_LIQ_aiNx", "UP_LIQ_aiNy", "UP_LIQ_aiNz", "VP_LIQ_aiNx", "VP_LIQ_aiNy", "VP_LIQ_aiNz", "WP_LIQ_aiNx", "WP_LIQ_aiNy", "WP_LIQ_aiNz"], 'P_extUaiN', 'z', r'$\langle Pu\delta^i \rangle$')
#pdu_l_ext=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dVdx","IP_LIQ_dWdx","IP_LIQ_dUdy", "IP_LIQ_dVdy","IP_LIQ_dWdy","IP_LIQ_dUdz","IP_LIQ_dVdz","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    pdu_l_ext_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dUdy","IP_LIQ_dUdz","IP_LIQ_dVdx", "IP_LIQ_dVdy","IP_LIQ_dVdz","IP_LIQ_dWdx","IP_LIQ_dWdy","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    p_l_ext_u_l_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UP_LIQ_I", "VP_LIQ_I", "WP_LIQ_I"],r'$Up_l$', r'$z$', r'$up_l$')

    vdvdxaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    duaiNx_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_i=Field.LoadFromFile(fstati,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    Idu_i = Field.LoadFromFile(fstati, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
####################################################################### Constructuion of the terms #################################################################################################
    ull_i=u_l_i*(I_i.inv(0.00001))
    uull_i = uu_l_i*(I_i.inv(0.00001))
 
    Rij_l_i=fluctuij(u_l_i, u_l_i, uu_l_i, I_i)
    uduain_i=vdvdxaiN_i+vdvdyaiN_i+vdvdzaiN_i  #Tensore che ha a che fare con il secondo termine viscoso all'interfaccia
    duain_i=dvdxaiN_i+dvdyaiN_i+dvdzaiN_i
    pression_ll_ext_i=pression_l_ext_i*(I_i.inv(0.00001)) 
    pu_int_l_i=p_l_ext_vaiN_i-pression_ll_ext_i*vaiN_c-ull_i.tensorProduct(p_l_ext_i)+pression_ll_ext_i*ull_i.tensorProduct(aiN_i)
    Idu_t_i = Idu_i.transpo()
    pdu_ll_ext_i=fluctuij(pression_l_ext_i, Idu_t_i, pdu_l_ext_i, I_i)
    pdu_ll_i = fluctuij(pression_l_i, Idu_t_i, pdu_l_i, I_i)
    pu_ll_ext_i=fluctuij(pression_l_ext_i, u_l_i, p_l_ext_u_l_i, I_i)
    pu_ll_i=fluctuij(pression_l_i, u_l_i, pu_l_i, I_i)
    dudu0_i=dudxdudx_l_i+dudydudy_l_i+dudzdudz_l_i ##This is a tensor

    dudu_l_i=fluctuij(Idu_i, Idu_i, dudu0_i, I_i) #Fluttuazione del gradiente di velocita 
    Rijk_l_i=fluctuijk(uuu_l_i, u_l_i, uu_l_i) #Prodotto delle velocita all'interno del termine di diffusione 

    Rij_i=Rij(u_l_i, ull_i, Rij_l_i, Rijk_l_i, pu_ll_ext_i, pdu_ll_ext_i, rho, dudu0_i, mu, I_i, uduain_i, duain_i, vaiN_i, aiN_i, vvaiN_i, pu_int_l_i, Idu_i)   
################################################################# second intermediate mesh ########################################################################################################
    dvar_2 = Field.getEntries(fstat2)
    I_2=Field.LoadFromFile(fstat2, ["coordonnee_K"] , ["I"],'I', 'z', r'I')
    u_l_2=Field.LoadFromFile(fstat2, ["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    uu_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    aiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    p_l_ext_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    pression_l_ext_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$P_l_ext$')
    uu_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    aii_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudxdudx_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdUdy", "IdUdxdUdx","IdUdxdUdy", "IdUdydUdy", "IdUdydUdz","IdUdxdUdz", "IdUdydUdz", "IdUdzdUdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudydudy_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdVdxdVdx", "IdVdxdVdy", "IdVdxdVdz","IdVdxdVdy", "IdVdydVdy", "IdVdydVdz","IdVdxdVdz", "IdVdydVdz", "IdVdzdVdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudzdudz_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdWdxdWdx", "IdWdxdWdy", "IdWdxdWdz","IdWdxdWdy", "IdWdydWdy", "IdWdydWdz","IdWdxdWdz", "IdWdydWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    p_l_ext_vaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"], ["UP_LIQ_aiNx", "UP_LIQ_aiNy", "UP_LIQ_aiNz", "VP_LIQ_aiNx", "VP_LIQ_aiNy", "VP_LIQ_aiNz", "WP_LIQ_aiNx", "WP_LIQ_aiNy", "WP_LIQ_aiNz"], 'P_extUaiN', 'z', r'$\langle Pu\delta^i \rangle$')
#pdu_l_ext=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dVdx","IP_LIQ_dWdx","IP_LIQ_dUdy", "IP_LIQ_dVdy","IP_LIQ_dWdy","IP_LIQ_dUdz","IP_LIQ_dVdz","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    pdu_l_ext_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dUdy","IP_LIQ_dUdz","IP_LIQ_dVdx", "IP_LIQ_dVdy","IP_LIQ_dVdz","IP_LIQ_dWdx","IP_LIQ_dWdy","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    p_l_ext_u_l_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UP_LIQ_I", "VP_LIQ_I", "WP_LIQ_I"],r'$Up_l$', r'$z$', r'$up_l$')

    vdvdxaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    duaiNx_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_2=Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    Idu_2 = Field.LoadFromFile(fstat2, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
####################################################################### Constructuion of the terms #################################################################################################
    ull_2=u_l_2*(I_2.inv(0.00001))
    uull_2 = uu_l_2*(I_2.inv(0.00001))
 
    Rij_l_2=fluctuij(u_l_2, u_l_2, uu_l_2, I_2)
    uduain_2=vdvdxaiN_2+vdvdyaiN_2+vdvdzaiN_2  #Tensore che ha a che fare con il secondo termine viscoso all'interfaccia
    duain_2=dvdxaiN_2+dvdyaiN_2+dvdzaiN_2
    pression_ll_ext_2=pression_l_ext_2*(I_2.inv(0.00001)) 
    pu_int_l_2=p_l_ext_vaiN_2-pression_ll_ext_2*vaiN_2-ull_2.tensorProduct(p_l_ext_2)+pression_ll_ext_2*ull_2.tensorProduct(aiN_2)
    Idu_t_2 = Idu_2.transpo()
    pdu_ll_ext_2=fluctuij(pression_l_ext_2, Idu_t_2, pdu_l_ext_2, I_2)
    pdu_ll_2 = fluctuij(pression_l_2, Idu_t_2, pdu_l_2, I_2)
    pu_ll_ext_2=fluctuij(pression_l_ext_2, u_l_2, p_l_ext_u_l_2, I_2)
    pu_ll_2=fluctuij(pression_l_2, u_l_2, pu_l_2, I_2)
    dudu0_2=dudxdudx_l_2+dudydudy_l_2+dudzdudz_l_2 ##This is a tensor

    dudu_l_2=fluctuij(Idu_2, Idu_2, dudu0_2, I_2) #Fluttuazione del gradiente di velocita 
    Rijk_l_2=fluctuijk(uuu_l_2, u_l_2, uu_l_2) #Prodotto delle velocita all'interno del termine di diffusione 

    Rij_2=Rij(u_l_2, ull_2, Rij_l_2, Rijk_l_2, pu_ll_ext_2, pdu_ll_ext_2, rho, dudu0_2, mu, I_2, uduain_2, duain_2, vaiN_2, aiN_2, vvaiN_2, pu_int_l_2, Idu_2)
####################################################################### Third mesh, indicated by f ###################################################################################################
    dvar_f = Field.getEntries(fstatf)
    I_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    u_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    uu_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    aiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["aiNx", "aiNy", "aiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    p_l_ext_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    pression_l_ext_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["PI"], 'PI', 'z' , r'$P_l_ext$')
    uu_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uuu_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    aii_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    pu_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UPI", "VPI", "WPI"],r'$Up_l$', r'$z$', r'$up_l$')
    dudxdudx_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdVdx", "IdUdxdWdx","IdUdxdVdx", "IdVdxdVdx", "IdVdxdWdx","IdUdxdWdx", "IdVdxdWdx", "IdWdxdWdx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudxdudx_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdUdxdUdx", "IdUdxdUdy", "IdUdxdUdx","IdUdxdUdy", "IdUdydUdy", "IdUdydUdz","IdUdxdUdz", "IdUdydUdz", "IdUdzdUdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudydudy_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["IdUdydUdy", "IdUdydVdy", "IdUdydWdy","IdUdydVdy", "IdVdydVdy", "IdVdydWdy","IdUdydWdy", "IdVdydWdy", "IdWdydWdy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudydudy_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdVdxdVdx", "IdVdxdVdy", "IdVdxdVdz","IdVdxdVdy", "IdVdydVdy", "IdVdydVdz","IdVdxdVdz", "IdVdydVdz", "IdVdzdVdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dudzdudz_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["IdUdzdUdz", "IdUdzdVdz", "IdUdzdWdz","IdUdzdVdz", "IdVdzdVdz", "IdVdzdWdz","IdUdzdWdz", "IdVdzdWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
#dudzdudz_l=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IdWdxdWdx", "IdWdxdWdy", "IdWdxdWdz","IdWdxdWdy", "IdWdydWdy", "IdWdydWdz","IdWdxdWdz", "IdWdydWdz", "IdWdzdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    pdu_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["IPdUdx","IPdUdy","IPdUdz","IPdVdx", "IPdVdy","IPdVdz","IPdWdx","IPdWdy","IPdWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UaiNx","UaiNy","UaiNz","VaiNx","VaiNy","VaiNz","WaiNx","WaiNy","WaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vPaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UPaiNx", "UPaiNy", "UPaiNz", "VPaiNx","VPaiNy", "VPaiNz", "WPaiNx", "WPaiNy", "WPaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vvaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUaiNx","UUaiNy","UUaiNz","UVaiNx","UVaiNy","UVaiNz","UWaiNx","UWaiNy","UWaiNz", "UVaiNx","UVaiNy","UVaiNz","VVaiNx","VVaiNy","VVaiNz","VWaiNx","VWaiNy","VWaiNz","UWaiNx","UWaiNy","UWaiNz","VWaiNx","VWaiNy","VWaiNz","WWaiNx","WWaiNy","WWaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    p_l_ext_vaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"], ["UP_LIQ_aiNx", "UP_LIQ_aiNy", "UP_LIQ_aiNz", "VP_LIQ_aiNx", "VP_LIQ_aiNy", "VP_LIQ_aiNz", "WP_LIQ_aiNx", "WP_LIQ_aiNy", "WP_LIQ_aiNz"], 'P_extUaiN', 'z', r'$\langle Pu\delta^i \rangle$')
#pdu_l_ext=Field.LoadFromFile(fstatc,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dVdx","IP_LIQ_dWdx","IP_LIQ_dUdy", "IP_LIQ_dVdy","IP_LIQ_dWdy","IP_LIQ_dUdz","IP_LIQ_dVdz","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    pdu_l_ext_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["IP_LIQ_dUdx","IP_LIQ_dUdy","IP_LIQ_dUdz","IP_LIQ_dVdx", "IP_LIQ_dVdy","IP_LIQ_dVdz","IP_LIQ_dWdx","IP_LIQ_dWdy","IP_LIQ_dWdz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$') 
    p_l_ext_u_l_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UP_LIQ_I", "VP_LIQ_I", "WP_LIQ_I"],r'$Up_l$', r'$z$', r'$up_l$')

    vdvdxaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UdUdxaiNx", "UdVdxaiNx", "UdWdxaiNx","VdUdxaiNx", "VdVdxaiNx", "VdWdxaiNx","WdUdxaiNx", "WdVdxaiNx", "WdWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdyaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UdUdyaiNy", "UdVdyaiNy", "UdWdyaiNy","VdUdyaiNy", "VdVdyaiNy", "VdWdyaiNy","WdUdyaiNy", "WdVdyaiNy", "WdWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    vdvdzaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UdUdzaiNz", "UdVdzaiNz", "UdWdzaiNz","VdUdzaiNz", "VdVdzaiNz", "VdWdzaiNz","WdUdzaiNz", "WdVdzaiNz", "WdWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdxaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    duaiNx_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_f=Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    Idu_f = Field.LoadFromFile(fstatf, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
###################################################################### Building variables ###############################################################################################
    ull_f=u_l_f*(I_f.inv(0.00001))
    uull_f = uu_l_f*(I_f.inv(0.00001))

    Rij_l_f=fluctuij(u_l_f, u_l_f, uu_l_f, I_f)
    uduain_f=vdvdxaiN_f+vdvdyaiN_f+vdvdzaiN_f  #Tensore che ha a che fare con il secondo termine viscoso all'interfaccia
    duain_f=dvdxaiN_f+dvdyaiN_f+dvdzaiN_f
    pression_ll_ext_f=pression_l_ext_f*(I_f.inv(0.00001)) 
    pu_int_l_f=p_l_ext_vaiN_f-pression_ll_ext_f*vaiN_f-ull_f.tensorProduct(p_l_ext_f)+pression_ll_ext_f*ull_f.tensorProduct(aiN_f)
    Idu_t_f = Idu_f.transpo()
    pdu_ll_ext_f=fluctuij(pression_l_ext_f, Idu_t_f, pdu_l_ext_f, I_f)
    pdu_ll_f = fluctuij(pression_l_f, Idu_t_f, pdu_l_f, I_f)
    pu_ll_ext_f=fluctuij(pression_l_ext_f, u_l_f, p_l_ext_u_l_f, I_f)
    pu_ll_f=fluctuij(pression_l_f, u_l_f, pu_l_f, I_f)
    dudu0_f=dudxdudx_l_f+dudydudy_l_f+dudzdudz_l_f ##This is a tensor

    dudu_l_f=fluctuij(Idu_f, Idu_f, dudu0_f, I_f) #Fluttuazione del gradiente di velocita 
    Rijk_l_f=fluctuijk(uuu_l_f, u_l_f, uu_l_f) #Prodotto delle velocita all'interno del termine di diffusione 
  
    Rij_f=Rij(u_l_f, ull_f, Rij_l_f, Rijk_l_f, pu_ll_ext_f, pdu_ll_ext_f, rho, dudu0_f, mu, I_f, uduain_f, duain_f, vaiN_f, aiN_f, vvaiN_f, pu_int_l_f, Idu_f)
    res ={'Rij_c':Rij_c,'Rij_i':Rij_i, 'Rij_f':Rij_f,'Rij_1':Rij_1, 'Rij_2':Rij_2}
    return res
######################################################################## Plotting function ####################################################
def Check(Residue2):
    Rij_c=Residue2['Rij_c']['interface']
    Rij_1=Residue2['Rij_1']['interface']
    Rij_i=Residue2['Rij_i']['interface']
    Rij_2=Residue2['Rij_2']['interface']
    Rij_f=Residue2['Rij_f']['interface']
    Rij_c.settex(r'1')
    Rij_i.settex(r'3')
    Rij_f.settex(r'5')
    tracer( [Rij_c, Rij_i, Rij_f], 'Res_x', [0], markevry=[1], markerSize=0.5)
    tracer( [Rij_c, Rij_i, Rij_f], 'Res_y', [4], markevry=[1], markerSize=0.5)
    tracer( [Rij_c, Rij_i, Rij_f], 'Res_z', [8], markevry=[1], markerSize=0.5)
    #a = Residue2['Rij_c']['pressure_int']
    p1 = Rij_c.compo(0,0)
    p1 = p1.int_mean()
    #b = Residue2['Rij_i']['pressure_int']
    p2 = Rij_i.compo(0,0)
    p2 = p2.int_mean()
    #c = Residue2['Rij_f']['pressure_int']
    p3 = Rij_f.compo(0,0)
    p3 = p3.int_mean()
    p4 = Rij_1.compo(0,0)
    p4 = p4.int_mean()
    p5 = Rij_2.compo(0,0)
    p5 = p5.int_mean()
    print p1, p2, p3, p4, p5 
    return


################################################################################### Function for Rij ##########################################################################################
def Rij(u_l, ull, Rij_l, Rijk_l, pu_ll_ext, pdu_ll_ext, rho, dudu_l, mu, I, uduain=None, duain=None, vaiN=None, aiN=None, vvaiN=None, pu_int_l=None, Idu=None):
    a5=u_l.tensorProduct(Rij_l)
    b5=a5.grad()
    inertia=b5.contract('ij','kijk')

    Idu_l = Idu*(I.inv(0.00001))
    nonlin2=Rij_l.tensorProduct(Idu_l)
    nonlin2=nonlin2.contract('ij','kijk')
    nonlin3=nonlin2.transpo()
    prod = (nonlin2+nonlin3) ####!!!!!!
    nonlin=inertia + prod 
################################ Turbulent diffusion ###############
    gradRijk=Rijk_l.grad()
    nonlin4=gradRijk.contract('ij','ikjk')
###################Pressure diffusion ##########
    pijk=pu_ll_ext.grad()    
    p1ijk=pijk.transpo()+pijk 
################# Redistibution #####################
    p2ijk=pdu_ll_ext.transpo()+pdu_ll_ext
################## Molecular diffusion #############
    gradRij=Rij_l.grad() 
    ddRij0=gradRij.grad()
    ddRij=ddRij0.contract('ij', 'ijkk')


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
    term2=term2.contract('ij', 'kijk') 
    term3=term3.contract('ij', 'kijk')

    viscinter2=(term1-term2+term3)  

    visc1=uduain + uduain.transpo() #Media totale valutata su tutto il dominio sfruttando la def di derivata del prodotto
    b=ull.tensorProduct(duain) 
    b= b + b.transpo()
###Attenzione cambiati termini di Idu ######
    c=Idu_l.tensorProduct(vaiN)
    c=c.contract('ij', 'kjik')
    c = c + c.transpo()
    visc2=b+c
 #termine legata al prodotto di u per la sua derivata di cui uno solo dei due risulta essere mediato 
    e=ull.tensorProduct(aiN)
    f=Idu_l.tensorProduct(e)    
    f_1=f.contract('ij', 'kjik')
    visc3 = f_1 + f_1.transpo()
    viscinter1=(visc1-visc2+visc3)


    advection=inertia*(-1.0)*(1./rho)
    production=(prod)*(-1.0)*(1./rho)
    redistribution=p2ijk*(1./rho)
    dissipation=dudu_l*mu*2.0*(-1.0)*(1./rho) 
    diffusion_turb=nonlin4*(-1.0) #+ ad1*(-1.0)
    transpression=p1ijk*(-1.0)*(1./rho) 
    diffusion_mol=ddRij*mu*(1./rho)
    diffusion = diffusion_turb + diffusion_mol + transpression
    interfacial_pressure=pu_int*(1./rho)
    viscinter = (viscinter1*(-1.0) +viscinter2*(-1.0))*mu*(1./rho)

    #for i in [dissipation, transpression, redistribution, diffusion_mol, interfacial_pressure, viscinter, diffusion_turb]:
    #        try:
    #            i=i*(1/rho)
    #        except:
    #            i=i*rho.inv()

    residu=(advection+production+redistribution+dissipation+transpression+diffusion_mol+diffusion_turb)*(-1.0) # e il RHS con cui confrontare il termine successivo
    m=interfacial_pressure+viscinter #segno presente all'interno dell'espressione 
    res_tot= residu - m

    
    res = {'advection':inertia, 'production': production, 'redistribution':redistribution, 'dissipation':dissipation, 'diffusion_mol':diffusion_mol, 'diffusion_turb':diffusion_turb, 'transpression':transpression, 'pressure_int':interfacial_pressure, 'viscous_int':viscinter, 'interface':m, 'RHS':residu, 'residual':res_tot}
    return res
Pippo=Residue2(fstatc, fstati, fstatf,fstat1, fstat2, rho, mu)
Check(Pippo)
