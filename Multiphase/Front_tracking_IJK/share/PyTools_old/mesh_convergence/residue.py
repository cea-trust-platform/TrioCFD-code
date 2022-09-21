# -*- coding: utf-8 -*-
#import rlcompleter
from Tools import *
import readline
import unittest
import pickle
import os


from matplotlib.pyplot import clf
readline.parse_and_bind('tab:complete') 
clf
global compteur
compteur=0
os.system('rm -f *.png')


gravity=-9.81  
rho_l=986.51
rho_v=845.44
alpha_v=0.0010472
mul=1.39e-3
muv=1.39e-3
Source=0.
sigma=0.00163  #sigma=0.00512    
fstatc = "Stats_coarse.dt_ev"   
timec = "Nsteady_coarse.dt_ev"
fstat1 = "Stats_int_1.dt_ev"   
time1 = "Nsteady_int_1.dt_ev"
fstat2 = "Stats_int_2.dt_ev"   
time2 = "Nsteady_int_2.dt_ev"
fstat3 = "Stats_int_3.dt_ev"
time3 = "Nsteady_3.dt_ev"
fstatf = "Stats_fin.dt_ev"
timef = "Nsteady_fin.dt_ev"



def Residue(fstatc,fstat1,fstat2,fstat3,fstatf,timec,time1,time2,time3,timef,gravity, rho_l, rho_v, alpha_v, mul, muv, Source, sigma):
    alpha_l=1.-alpha_v
    rho_ll=rho_l*alpha_l
    rho_vv=rho_v*alpha_v
    rho_av=rho_ll+rho_vv
    
########################################################## pression_interfacialiable_coarse_mesh (identified by a letter c) #################################################################################### 
    dpression_interfacialc = Field.getEntries(fstatc)
    I_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    pression_l_ext_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_v_ext_c = Field.LoadFromFile(fstatc,["coordonnee_K"], ["P_VAP_Iv"], 'P_VAP_Iv', 'z', r'$P_v_ext$')
    p_l_ext_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    p_v_ext_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],'p_v', r'$z$', r'$p_v\nabla\chi_v$')#*(-1.0)
    uu_l_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    duaiNx_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    dvdxaiN_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_c = Field.LoadFromFile(fstatc,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    Idu_c = Field.LoadFromFile(fstatc, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
    Ivdu_c = Field.LoadFromFile(fstatc, ["coordonnee_K"], ["IvdUdx", "IvdVdx", "IvdWdx", "IvdUdy", "IvdVdy", "IvdWdy", "IvdUdz", "IvdVdz", "IvdWdz"], 'Ivdu', 'z', r'$\overline{\chi \nabla u}$')   

    Iv_c = (I_c-1.0)*(-1.0)

    g_c=Field.initgravity([gravity,0,0], 'g', I_c)
    S_c=Field.initsource([Source,0,0], 'beta', I_c)

    ########### Time derivative ############################# 
   
    dpression_interfacialcc = Field.getEntries(timec)
    u_l_cc=Field.LoadFromFile(timec,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    u_v_cc=Field.LoadFromFile(timec,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_v$', r'$z$', r'$u_v$') 
    ############################################################### Terms added for the momentum equation ###########################################################  
    p_v_ext_c=p_v_ext_c*(-1.0)   
    Rij_l_c = fluctuij(u_l_c, u_l_c, uu_l_c, I_c)
    Rij_v_c = fluctuij(u_v_c, u_v_c, uu_v_c, Iv_c)
    tau_liq_c = duaiNx_c + dvdxaiN_c + dvaiNy_c + dvdyaiN_c + dwaiNz_c + dvdzaiN_c
    tau_vap_c = tau_liq_c*(-1.0) 
    SRT_liq_c=Idu_c
    SRT_vap_c=Ivdu_c
    
    qdm_liq_c=qdm('liquide',u_l_c, mul, pression_l_ext_c,   Rij_l_c, S_c, rho_l, u_l_cc, aii_c, rho_av, g_c, sigma, I_c, p_l_ext_c, tau_liq_c, SRT_liq_c)
    qdm_gaz_c=qdm('gaz',u_v_c, muv, pression_v_ext_c,  Rij_v_c, S_c, rho_v, u_v_cc, aii_c, rho_av, g_c, sigma, Iv_c, p_v_ext_c, tau_vap_c, SRT_vap_c) 

#########################################################pression_interfacialiable inter mesh 1(identified by the number 1) ######################################################################
    dpression_interfacial1 = Field.getEntries(fstat1)
    I_1 =Field.LoadFromFile(fstat1,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    pression_l_ext_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_v_ext_1 = Field.LoadFromFile(fstat1,["coordonnee_K"], ["P_VAP_Iv"], 'P_VAP_Iv', 'z', r'$P_v_ext$')
    p_l_ext_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    p_v_ext_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],'p_v', r'$z$', r'$p_v\nabla\chi_v$')#*(-1.0)
    uu_l_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    duaiNx_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    dvdxaiN_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_1 = Field.LoadFromFile(fstat1,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    Idu_1 = Field.LoadFromFile(fstat1, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
    Ivdu_1 = Field.LoadFromFile(fstat1, ["coordonnee_K"], ["IvdUdx", "IvdVdx", "IvdWdx", "IvdUdy", "IvdVdy", "IvdWdy", "IvdUdz", "IvdVdz", "IvdWdz"], 'Ivdu', 'z', r'$\overline{\chi \nabla u}$')

    Iv_1 = (I_1-1.0)*(-1.0)

    g_1=Field.initgravity([gravity,0,0], 'g', I_1)
    S_1=Field.initsource([Source,0,0], 'beta', I_1)

    dpression_interfacial11 = Field.getEntries(time1)
    u_l_11=Field.LoadFromFile(time1,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    u_v_11=Field.LoadFromFile(time1,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_v$', r'$z$', r'$u_v$') 
    
    
      
    p_v_ext_1=p_v_ext_1*(-1.0)   
    Rij_l_1 = fluctuij(u_l_1, u_l_1, uu_l_1, I_1)
    Rij_v_1 = fluctuij(u_v_1, u_v_1, uu_v_1, Iv_1)
    tau_liq_1 = duaiNx_1 + dvdxaiN_1 + dvaiNy_1 + dvdyaiN_1 + dwaiNz_1 + dvdzaiN_1
    tau_vap_1 = tau_liq_1*(-1.0) 
    SRT_liq_1=Idu_1
    SRT_vap_1=Ivdu_1
    

    qdm_liq_1=qdm('liquide',u_l_1, mul, pression_l_ext_1,   Rij_l_1, S_1, rho_l, u_l_11, aii_1, rho_av, g_1, sigma, I_1, p_l_ext_1, tau_liq_1, SRT_liq_1)
    qdm_gaz_1=qdm('gaz',u_v_1, muv, pression_v_ext_1,  Rij_v_1, S_1, rho_v, u_v_11, aii_1, rho_av, g_1, sigma, Iv_1, p_v_ext_1, tau_vap_1, SRT_vap_1) 
  
    
#########################################################pression_interfacialiable inter mesh 2(identified by the number 2) ######################################################################
    dpression_interfacial2 = Field.getEntries(fstat2)
    I_i =Field.LoadFromFile(fstat2,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    pression_l_ext_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_v_ext_i = Field.LoadFromFile(fstat2,["coordonnee_K"], ["P_VAP_Iv"], 'P_VAP_Iv', 'z', r'$P_v_ext$')
    p_l_ext_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    p_v_ext_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],'p_v', r'$z$', r'$p_v\nabla\chi_v$')#*(-1.0)
    uu_l_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    duaiNx_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    dvdxaiN_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_i = Field.LoadFromFile(fstat2,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    Idu_i = Field.LoadFromFile(fstat2, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
    Ivdu_i = Field.LoadFromFile(fstat2, ["coordonnee_K"], ["IvdUdx", "IvdVdx", "IvdWdx", "IvdUdy", "IvdVdy", "IvdWdy", "IvdUdz", "IvdVdz", "IvdWdz"], 'Ivdu', 'z', r'$\overline{\chi \nabla u}$')

    Iv_i = (I_i-1.0)*(-1.0)

    g_i=Field.initgravity([gravity,0,0], 'g', I_i)
    S_i=Field.initsource([Source,0,0], 'beta', I_i)

    dpression_interfacial22 = Field.getEntries(time2)
    u_l_2=Field.LoadFromFile(time2,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    u_v_2=Field.LoadFromFile(time2,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_v$', r'$z$', r'$u_v$') 
    #u_2=u_l_2*rho_l+u_v_2*rho_v
    
    #u_i=u_v_i+u_l_i #velocity of the mixture coarse mesh#   
    p_v_ext_i=p_v_ext_i*(-1.0)   
    #uu=uu_v+uu_l
    #pression=pression_l+pression_v 
    #pression_new=pression_l_ext+pression_v_ext        
    #rhog_av=g*rho_av
    #rhog_tmp=rho*g
    #rhog=rhog_tmp-rhog_av
    #ufu=u.MatProduct(u)
    #Rij=uu-ufu           #media del prodotto meno prodotto delle medie 
    #ufu_l=u_l.MatProduct(u_l)
    #Rij_l=uu_l-ufu_l
    Rij_l_i = fluctuij(u_l_i, u_l_i, uu_l_i, I_i)
    #ufu_v=u_v.MatProduct(u_v)
    #Rij_v=uu_v-ufu_v
    Rij_v_i = fluctuij(u_v_i, u_v_i, uu_v_i, Iv_i)
    tau_liq_i = duaiNx_i + dvdxaiN_i + dvaiNy_i + dvdyaiN_i + dwaiNz_i + dvdzaiN_i
    tau_vap_i = tau_liq_i*(-1.0) 
    SRT_liq_i=Idu_i
    SRT_vap_i=Ivdu_i
    

    qdm_liq_i=qdm('liquide',u_l_i, mul, pression_l_ext_i,   Rij_l_i, S_i, rho_l, u_l_2, aii_i, rho_av, g_i, sigma, I_i, p_l_ext_i, tau_liq_i, SRT_liq_i)
    qdm_gaz_i=qdm('gaz',u_v_i, muv, pression_v_ext_i,  Rij_v_i, S_i, rho_v, u_v_2, aii_i, rho_av, g_i, sigma, Iv_i, p_v_ext_i, tau_vap_i, SRT_vap_i) 
 

#########################################################pression_interfacialiable inter 3 mesh (identified by the number 3) ######################################################################
    dpression_interfacial3 = Field.getEntries(fstat3)
    I_3 =Field.LoadFromFile(fstat3,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    pression_l_ext_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_v_ext_3 = Field.LoadFromFile(fstat3,["coordonnee_K"], ["P_VAP_Iv"], 'P_VAP_Iv', 'z', r'$P_v_ext$')
    p_l_ext_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    p_v_ext_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],'p_v', r'$z$', r'$p_v\nabla\chi_v$')#*(-1.0)
    uu_l_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    duaiNx_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    dvdxaiN_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_3 = Field.LoadFromFile(fstat3,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    Idu_3 = Field.LoadFromFile(fstat3, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
    Ivdu_3 = Field.LoadFromFile(fstat3, ["coordonnee_K"], ["IvdUdx", "IvdVdx", "IvdWdx", "IvdUdy", "IvdVdy", "IvdWdy", "IvdUdz", "IvdVdz", "IvdWdz"], 'Ivdu', 'z', r'$\overline{\chi \nabla u}$')

    Iv_3 = (I_3-1.0)*(-1.0)

    g_3=Field.initgravity([gravity,0,0], 'g', I_3)
    S_3=Field.initsource([Source,0,0], 'beta', I_3)


    dpression_interfacial33 = Field.getEntries(time3)
    u_l_33=Field.LoadFromFile(time3,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    u_v_33=Field.LoadFromFile(time3,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_v$', r'$z$', r'$u_v$') 
    #u_2=u_l_2*rho_l+u_v_2*rho_v
    
    #u_i=u_v_i+u_l_i #velocity of the mixture coarse mesh#   
    p_v_ext_3=p_v_ext_3*(-1.0)   
    #uu=uu_v+uu_l
    #pression=pression_l+pression_v 
    #pression_new=pression_l_ext+pression_v_ext        
    #rhog_av=g*rho_av
    #rhog_tmp=rho*g
    #rhog=rhog_tmp-rhog_av
    #ufu=u.MatProduct(u)
    #Rij=uu-ufu           #media del prodotto meno prodotto delle medie 
    #ufu_l=u_l.MatProduct(u_l)
    #Rij_l=uu_l-ufu_l
    Rij_l_3 = fluctuij(u_l_3, u_l_3, uu_l_3, I_3)
    #ufu_v=u_v.MatProduct(u_v)
    #Rij_v=uu_v-ufu_v
    Rij_v_3 = fluctuij(u_v_3, u_v_3, uu_v_3, Iv_3)
    tau_liq_3 = duaiNx_3 + dvdxaiN_3 + dvaiNy_3 + dvdyaiN_3 + dwaiNz_3 + dvdzaiN_3
    tau_vap_3 = tau_liq_3*(-1.0) 
    SRT_liq_3=Idu_3
    SRT_vap_3=Ivdu_3
    

    qdm_liq_3=qdm('liquide',u_l_3, mul, pression_l_ext_3,   Rij_l_3, S_3, rho_l, u_l_33, aii_3, rho_av, g_3, sigma, I_3, p_l_ext_3, tau_liq_3, SRT_liq_3)
    qdm_gaz_3=qdm('gaz',u_v_3, muv, pression_v_ext_3,  Rij_v_3, S_3, rho_v, u_v_33, aii_3, rho_av, g_3, sigma, Iv_3, p_v_ext_3, tau_vap_3, SRT_vap_3) 

    
#########################################################pression_interfacialiable final mesh (identified by the letter f) ######################################################################
    dpression_interfacialf = Field.getEntries(fstatf)
    I_f =Field.LoadFromFile(fstatf,["coordonnee_K"] , ["I"],'I', 'z', r'I')
    pression_l_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["PI"],'P_l', 'z', r'$P_l$')
    pression_v_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["PIv"],'P_v', 'z', r'$P_v$')
    pression_l_ext_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["P_LIQ_I"], 'P_LIQ_I', 'z' , r'$P_l_ext$')
    pression_v_ext_f = Field.LoadFromFile(fstatf,["coordonnee_K"], ["P_VAP_Iv"], 'P_VAP_Iv', 'z', r'$P_v_ext$')
    p_l_ext_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],'p_l', r'$z$', r'$p_l\nabla\chi_l$')
    p_v_ext_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],'p_v', r'$z$', r'$p_v\nabla\chi_v$')#*(-1.0)
    uu_l_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUI", "UVI", "UWI","UVI", "VVI", "VWI","UWI", "VWI", "WWI"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    uu_v_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUIv", "UVIv", "UWIv","UVIv", "VVIv", "VWIv","UWIv", "VWIv", "WWIv"],'UU_v', 'z', r'\langle(u\times u)\rangle)_v')
    uuu_l_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UUUI", "UUVI", "UUWI","UUVI", "UVVI", "UVWI","UUWI", "UVWI", "UWWI","UUVI", "UVVI", "UVWI","UVVI", "VVVI", "VVWI","UVWI", "VVWI", "VWWI","UUWI", "UVWI", "UWWI","UVWI", "VVWI", "VWWI","UWWI", "VWWI", "WWWI" ],'UUU', 'z', r'\langle(u\times u\times u)\rangle)_v')
    u_l_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$')
    u_v_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["UIv", "VIv", "WIv"],'U_v', 'z', r'$u_v$')
    aii_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["kaiNx", "kaiNy", "kaiNz"],'f_interf', 'z', r'$\sigma n\kappa\delta^i$')
    duaiNx_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdxaiNx","dUdyaiNx","dUdzaiNx"], 'UU_l', 'z', r'$\nabla\cdot\u_l$')
    dvaiNy_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dVdxaiNy","dVdyaiNy","dVdzaiNy"], 'UU_l', 'z', r'$\nabla\cdot\v_l$')
    dwaiNz_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dWdxaiNz","dWdyaiNz","dWdzaiNz"], 'UU_l', 'z', r'$\nabla\cdot\w_l$')
    dvdxaiN_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdyaiN_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    dvdzaiN_f = Field.LoadFromFile(fstatf,["coordonnee_K"] , ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],'UU_l', 'z', r'$\langle(u\times u)\rangle)_l$')
    Idu_f = Field.LoadFromFile(fstatf, ["coordonnee_K"], ["IdUdx", "IdVdx", "IdWdx", "IdUdy", "IdVdy", "IdWdy", "IdUdz", "IdVdz", "IdWdz"], 'Idu', 'z', r'$\overline{\chi \nabla u}$')
    Ivdu_f = Field.LoadFromFile(fstatf, ["coordonnee_K"], ["IvdUdx", "IvdVdx", "IvdWdx", "IvdUdy", "IvdVdy", "IvdWdy", "IvdUdz", "IvdVdz", "IvdWdz"], 'Ivdu', 'z', r'$\overline{\chi \nabla u}$')

    Iv_f = (I_f-1.0)*(-1.0)

    g_f=Field.initgravity([gravity,0,0], 'g', I_f)
    S_f=Field.initsource([Source,0,0], 'beta', I_f)


    dpression_interfacialff = Field.getEntries(timef)
    u_l_ff=Field.LoadFromFile(timef,["coordonnee_K"] , ["UI", "VI", "WI"],r'$U_l$', r'$z$', r'$u_l$') 
    u_v_ff=Field.LoadFromFile(timef,["coordonnee_K"] , ["UIv", "VIv", "WIv"],r'$U_v$', r'$z$', r'$u_v$') 
    #u_2=u_l_2*rho_l+u_v_2*rho_v
    
    #u_i=u_v_i+u_l_i #velocity of the mixture coarse mesh#   
    p_v_ext_f=p_v_ext_f*(-1.0)
    Rij_l_f = fluctuij(u_l_f, u_l_f, uu_l_f, I_f)
    Rij_v_f = fluctuij(u_v_f, u_v_f, uu_v_f, Iv_f)
    tau_liq_f = duaiNx_f + dvdxaiN_f + dvaiNy_f + dvdyaiN_f + dwaiNz_f + dvdzaiN_f
    tau_vap_f = tau_liq_f*(-1.0) 
    SRT_liq_f=Idu_f
    SRT_vap_f=Ivdu_f
    

    qdm_liq_f=qdm('liquide',u_l_f, mul, pression_l_ext_f,   Rij_l_f, S_f, rho_l, u_l_ff, aii_f, rho_av, g_f, sigma, I_f, p_l_ext_f, tau_liq_f, SRT_liq_f)
    qdm_gaz_f=qdm('gaz',u_v_f, muv, pression_v_ext_f,  Rij_v_f, S_f, rho_v, u_v_ff, aii_f, rho_av, g_f, sigma, Iv_f, p_v_ext_f, tau_vap_f, SRT_vap_f) 
    res={'qdm_liq_i':qdm_liq_i,'qdm_liq_c':qdm_liq_c,'qdm_gaz_c':qdm_gaz_c,'qdm_gaz_i':qdm_gaz_i, 'qdm_liq_3':qdm_liq_3, 'qdm_gaz_3':qdm_gaz_3, 'qdm_liq_1':qdm_liq_1, 'qdm_gaz_1':qdm_gaz_1, 'qdm_liq_f':qdm_liq_f, 'qdm_gaz_f':qdm_gaz_f}

    return res

   ############################################################# Plot of the Interfaceious quantities ############################################################# 
def Check(Residue):  #### Change the argument to have the confront between a different quantity ###################

    Aliq_c=Residue['qdm_liq_c']['Interface']
    Avap_c=Residue['qdm_gaz_c']['Interface']
    
    #Aliq_c=Aliq_c/mul
    #Avap_c=Avap_c/muv

    Aliq_1=Residue['qdm_liq_1']['Interface']
    Avap_1=Residue['qdm_gaz_1']['Interface']
    
    #Aliq_1=Aliq_1/mul
    #Avap_1=Avap_1/muv

    Aliq_i=Residue['qdm_liq_i']['Interface']
    Avap_i=Residue['qdm_gaz_i']['Interface']
    
    #Aliq_i=Aliq_i/mul
    #Avap_i=Avap_i/muv
    
    Aliq_3=Residue['qdm_liq_3']['Interface']
    Avap_3=Residue['qdm_gaz_3']['Interface']
    
    #Aliq_3=Aliq_3/mul
    #Avap_3=Avap_3/muv
    
    Aliq_f=Residue['qdm_liq_f']['Interface']
    Avap_f=Residue['qdm_gaz_f']['Interface']
    
    #Aliq_f=Aliq_f/mul
    #Avap_f=Avap_f/muv


####################### To do only on the refined mesh to check the balance of the two phases ########################################
####################### For integrated field ##################################
    #M_sigma=Residue['qdm_liq_f']['Sigma']*(-1.0)
    #Ml=Residue['qdm_liq_f']['M']
    #Mv=Residue['qdm_gaz_f']['M']
####################### For the instantaneous field ######################################
    #M_sigma=Residue['qdm_liq_f']['sigma']*(-1.0)
    #Ml=Residue['qdm_liq_f']['Interface']
    #Mv=Residue['qdm_gaz_f']['Interface']

####################### For the slides #######################################
    M_sigma=Residue['qdm_liq_f']['sigma']#*(-1.0)
    Mlp=Residue['qdm_liq_f']['pression_interfacial']
    Mlv=Residue['qdm_liq_f']['Stress_tensor']
    Mvp=Residue['qdm_gaz_f']['pression_interfacial']
    Mvv=Residue['qdm_gaz_f']['Stress_tensor']
    M_sigma_post = Mlp + Mlv + Mvp + Mvv 
    ########## Print settings ############################################
    #M_sigma_post=Mv+Ml
    ########## Print settings ############################################
    #Ml.settex(r'$\int_{0}^{y} \mathbf{M_{l}} dy$ ')
    #Mv.settex(r'$\int_{0}^{y} \mathbf{M_{v}} dy$')
    #M_sigma.settex(r'$\int_{0}^{y} \mathbf{M_{\sigma}} dy$')
    M_sigma.settex(r'$\mathbf{M_{\sigma}}$')
    #M_sigma_post.settex(r'$\int_{0}^{y} \mathbf{M_{l}} + \mathbf{M_{v}} dy$')
    M_sigma_post.settex(r'$\mathbf{M_{l}} + \mathbf{M_{v}}$')
    Mlp.settex(r'$\mathbf{M^{p}_{l}}$')
    Mvp.settex(r'$\mathbf{M^{p}_{v}}$')
    Mlv.settex(r'$\mathbf{M^{\mu}_{l}}$')
    Mvv.settex(r'$\mathbf{M^{\mu}_{v}}$')
    tracer([Mlp, Mlv, Mvp, Mvv, M_sigma, M_sigma_post], 'Mx', [0], markevry=[0.04])
    tracer([Mlp, Mlv, Mvp, Mvv, M_sigma, M_sigma_post], 'Mz', [2], markevry=[0.04], markerSize=0.5)
    #M_sigma_old.settex(r'$M_{l}+M_{v}$')
    #Mll.settex(r'$M_{ll}$')
    #Mvv.settex(r'$M_{vv}$')
    #tracer([Mll, Mvv, M_sigma, M_sigma_old], 'Mx', [0], markevry=[1], markerSize=0.5)
    #tracer([Mll, Mvv, M_sigma, M_sigma_old], 'My', [1], markevry=[1], markerSize=0.5 )
    #tracer([Mll, Mvv, M_sigma, M_sigma_old], 'Mz', [2], markevry=[1], markerSize=0.5)
    #tracer([Ml, Mv, M_sigma,  M_sigma_post], 'Mx_new', [0], markevry=[0.025])
    #tracer([Ml, Mv, M_sigma,  M_sigma_post], 'My_new', [1], markevry=[0.025])
    #tracer([Ml, Mv, M_sigma,  M_sigma_post], 'Mz_new', [2], markevry=[0.025])

######################################################################################################################################################################## 
    
    Aliq_c.settex(r'M_1_liquide')
    Aliq_1.settex(r'M_2_liquide')
    Aliq_i.settex(r'M_3_liquide')
    Aliq_3.settex(r'M_4_liquide')
    Aliq_f.settex(r'M_5_liquide')
    
    Avap_c.settex(r'M_1_gas')
    Avap_1.settex(r'M_2_gas') 
    Avap_i.settex(r'M_3_gas')
    Avap_3.settex(r'M_4_gas')
    Avap_f.settex(r'M_5_gas')

    tracer([Aliq_c, Aliq_1, Aliq_i, Aliq_3, Aliq_f], 'Res_liq_x', [0], markevry=[1], markerSize=0.5)
    tracer([Aliq_c, Aliq_1,  Aliq_i, Aliq_3, Aliq_f],'Res_liq_y', [1], markevry=[1], markerSize=0.5)
    tracer([Aliq_c, Aliq_1, Aliq_i, Aliq_3, Aliq_f], 'Res_liq_z', [2], markevry=[1], markerSize=0.5)
    
    tracer([Avap_c, Avap_1, Avap_i, Avap_3, Avap_f], 'Res_vap_x', [0], markevry=[1], markerSize=0.5)
    tracer([Avap_c, Avap_1, Avap_i, Avap_3, Avap_f], 'Res_vap_y', [1], markevry=[1], markerSize=0.5)
    tracer([Avap_c, Avap_1, Avap_i, Avap_3, Avap_f], 'Res_vap_z', [2], markevry=[1], markerSize=0.5)


    #e_1 = 48./0.005
    #e_2 = 72./0.005
    #e_3 = 96./0.005
    #e_4 = 136/0.005
    #e_5 = 192/0.005
        
    ####################### Provides the mean values: use int_mean for interfacial quantity only, b,c,d if you want RMS
    a_1=Aliq_c.compo(0,0)
    #b_1 = a_1.power(2)
    #c_1 = b_1.mean()
    #d_1 = np.sqrt(c_1)
    a_1 = a_1.int_mean()
    a_2=Aliq_1.compo(0,0)
    #b_2 = a_2.power(2)
    #c_2 = b_2.mean()
    #d_2 = np.sqrt(c_2)
    a_2 = a_2.int_mean()
    a_3=Aliq_i.compo(0,0)
    #b_3 = a_3.power(2)
    #c_3 = b_3.mean()
    #d_3 = np.sqrt(c_3)
    a_3=a_3.int_mean()
    a_4=Aliq_3.compo(0,0)
    a_4= a_4.int_mean()
    #b_4 = a_4.power(2)
    #c_4 = b_4.mean()
    #d_4 = np.sqrt(c_4)
    a_5=Aliq_f.compo(0,0)
    a_5= a_5.int_mean()
    #b_5 = a_5.power(2)
    #c_5 = b_5.mean()
    #d_5 = np.sqrt(c_5)
    
    
     
    print d_1, d_2, d_3, d_4, d_5
    return

    ##################################################### Momentum equation ##############################################################################################
     
def qdm(phase, u, mu, pression, Rij, S, rho, u_1, aii=None, rho_av=None, g=None, sigma=None, I=None, pext=None, tau_int=None, SRT=None):
    d = aii*sigma
    D = d.integ()
    tau1=SRT
    tau2=tau1+tau1.transpo()
    tau3=tau2.grad()
    tau4=tau3.contract('i', 'ikk') #this is the divergence of the shear stress tensor  
    nonlin=u.tensorProduct(u)
    gradp=pression.grad() #The field is the average in the one-fluid equation
    divuu=nonlin.grad()
    divuu=divuu.contract('i', 'ikk') #Divergence of the convective term fro the averaged field. Does this term appear at the statisctical equilibrium???
    divRij=Rij.grad()
    divRij=divRij.contract('i', 'ikk')
    nsteady = u_1*rho*(-1.0) 
    divuu=divuu*rho
    divRij=divRij*rho
    divT=tau4*mu      
    tau_int=tau_int*mu
    m=pext*(-1.0)+tau_int #Segno corretto se il termine e preso come right hand side
    rhog_av=g*rho_av
    rhog_tmp=g*rho
    rhog=rhog_tmp-rhog_av
    rhog=I*rhog
    S=I*S
    ai=divT+divRij*(-1.0)+gradp*(-1.0)+S+divuu*(-1.0)+rhog + nsteady  
    gradp=gradp*(-1.0)
    divRij=divRij*(-1.0)
    divuu=divuu*(-1.0)
######## For the last plot only
    P_ext=pext.integ()
    Tau_ext=tau_int.integ()
##########################3 
    tau_int = tau_int*(-1.0)
    insta=ai+pext+tau_int  #new total residue of the equation 
    M =m.integ()
    P_ext=pext.integ()
    Tau_ext=tau_int.integ()
    res={'inertie':divuu, 'turbulence':divRij, 'viscous':divT, 'pression':gradp, 'source':S, 'RHS':ai, 'rhog':rhog, 'pression_interfacial':pext, 'Stress_tensor':tau_int, 'Interface':m, 'residue':insta, 'Sigma':D, 'M':M, 'sigma':d, 'Rij': Rij, 'P_ext_int':P_ext, 'Tau_ext_int':Tau_ext}
    return res

Pippo=Residue(fstatc,fstat1,fstat2,fstat3,fstatf,timec,time1,time2,time3,timef,gravity, rho_l, rho_v, alpha_v, mul, muv, Source, sigma)
Check(Pippo)
