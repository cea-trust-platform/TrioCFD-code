# -*- coding: utf-8 -*-
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import glob
import sys
import DNSTools3 as dtool

################################################################################################
# -> A bouger dans DNSTools je pense...
# How to correctly sort a string with a number inside? [duplicate] ;
# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
import re

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]
################################################################################################


############################################################## Mode d'emploi #######################################
# mode d'emploi :  python vitesses_phases.py who_to_play_with nom_des_figures liste_des_fichier.data liste_des_chemins_vers_diphasique_moyenne_spatiale**.txt liste_des_signes [t_deb,t_fin] corr
#                                 0              1                2                 3                                    4                                           5              6         7
# /!\ 1 : who_to_play_with = xyz_liq-vap -> trace les composantes x,y et z de la vitesse pour la phase liquie et vapeur le tout sur le meme graphique
# /!\ 1 : who_to_play_with = x_vap -> trace la composante x de la vitesse pour la phase liquie le tout sur le meme graphique
# /!\ 1 : who_to_play_with = all_all -> trace toutes les composantes de la vitesse pour toutes les phases, le tout sur le meme graphique
# /!\ 3 : [CHEMIN_data1, CHEMIN_data2, ...]les chemins DOIVENT Se terminer par .data
# /!\ 4 : [CHEMIN_TXT1, CHEMIN_TXT2, ...]les chemins DOIVENT Se terminer par un "/"
# /!\ 5 : [temps0,temps1] -> releve les donnees des txt entre les 2 temps
# /!\ 5 : [temps0,temps0] -> releve la donnee des outs a l'instant le plus proche de temps0
#         [-1,-1] -> releve tous les temps disponibles
# 7 : mets 1 si tu as mis BIG et 10 pour post_trt_derive.sh. Met 0 autrement
# Exemple sur gr262753@is243146:
# /volatile/RANDOM_TRIPERIO$ python /volatile/RANDOM_TRIPERIO/SOME_PYTHON/vitesses_phases_fusion_cas.py yz_liq fig_U_mix [petit_a015/COARSE/MERGE_0_4/DNS_swarm.data,petit_a030/COARSE/MERGE_0_4/DNS_swarm.data,petit_a060/COARSE/RUN01/DNS_swarm.data,petit_a120/COARSE/MERGE_0_6/DNS_swarm.data] [petit_a015/COARSE/MERGE_0_4/,petit_a030/COARSE/MERGE_0_4/,petit_a060/COARSE/MERGE_0_6/,petit_a120/COARSE/MERGE_0_6/] [-1,-1]
############################################################## fin Mode d'emploi #######################################
# Remarque generale :
# L'affichage des pentes, et des courbes "corrigee" n'est proprement fait QUE POUR la phase liquide !!
# Il y a un début pour la phase vapeur, mais c'est a relire entierement !!
############################################################## Mode d'emploi #######################################
liste_chemin = sys.argv[4].rstrip("]").lstrip("[").split(",")
liste_jdd = sys.argv[3].rstrip("]").lstrip("[").split(",")
good_liste_txt = []
good_liste_data = []
# # Rustine a la main parce qu'on veut mettre bout a bout des N00 et des N01 pour tous les taux de vide. Sauf le 12% qui n'a pas de N01...
# BIG_dix = int(sys.argv[7])
# print("BIG_dix",BIG_dix)
# if BIG_dix:
#   for txt in liste_chemin:
#     if (("N01" in txt) or ("SIMU_120" in txt)) :
#         good_liste_txt.append(txt)
#   for data in liste_jdd:
#     if (("N01" in data) or ("SIMU_120" in data)) : 
#         good_liste_data.append(data)
# 
#   liste_chemin = good_liste_txt
#   liste_jdd = good_liste_data



liste_signes = sys.argv[5].rstrip("]").lstrip("[").split(",")
print("#######################")
print(liste_chemin)
print(liste_jdd)
# fstat = sys.argv[4]
composantes, phases = sys.argv[1].split("_")
t_deb, t_fin = float(sys.argv[6].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[6].rstrip("]").lstrip("[").split(",")[-1])
figname = sys.argv[2]


print(composantes)
print(t_deb,t_fin)

if (t_deb <= 0.):
    t_deb = 0
if (t_fin <= 0.):
    t_fin = 1e6

# On ouvre la figure
fig, ax = plt.subplots(1,figsize=(14,7))
fig_smt, ax_smt = plt.subplots(1,figsize=(14,7))
fig_corr, ax_corr = plt.subplots(1,figsize=(14,7))
fig_cp, ax_cp = plt.subplots(1,figsize=(14,7))
fig_pente, ax_pente = plt.subplots(1,figsize=(14,7))

ax.set_title('Vitesses, moyenne de phase '+r' : $\overline{U_i}^k = \frac{\overline{\chi_k U_i}}{\overline{\chi_k}}$',fontsize=18)
ax_smt.set_title('Vitesses, moyenne de phase lissees en temps'+r' : $\overline{U_i}^k = \frac{\overline{\chi_k U_i}}{\overline{\chi_k}}$',fontsize=18)
ax_corr.set_title('Vitesses, moyenne de phase '+r' : $\overline{U_i}^k = \frac{\overline{\chi_k U_i}}{\overline{\chi_k}}$, corrigees',fontsize=18)
ax_cp.set_title('Vitesses, moyenne de phase '+r' : $\overline{U_i}^k = \frac{\overline{\chi_k U_i}}{\overline{\chi_k}}$, cor+pentes',fontsize=18)
ax_pente.set_title('Vitesses, moyenne de phase '+r' : $\overline{U_i}^k = \frac{\overline{\chi_k U_i}}{\overline{\chi_k}}$, pentes',fontsize=18)

linestyle = ["solid","dashed","dashdot","dotted"]
if len(liste_chemin)>0:
    linestyle_rev = linestyle 
    linestyle_rev.reverse()
linewidth = [1.,2.,3.,4.]
alpha = [1.0,0.6,0.4,0.25]
# linestyle2 : loosely dotted; loosely dashed; loosely dashdotted; solid (ofc
linestyle2 = [(0, (1,10)), (0, (5, 10)), (0, (3,10,1,10)), "solid"]
if len(liste_chemin)==1:
    linestyle2.reverse()
for ind, jdd in enumerate(liste_jdd) :
    # fstat = liste_chemin[ind]
    
    mu = dtool.getParam(jdd, "mu_liquide")
    rho = dtool.getParam(jdd, "rho_liquide")
    nu = mu/rho
    print("######################################")
    print("# "+str(liste_signes[ind]))
    print("## parametres who_to_play_with : "+composantes+", "+phases)
    if ("all" in composantes):
        composantes = "xyz"
    if ("all" in phases):
        phases = "liq-vap"
    print("# composantes de vitesses : "+composantes)
    print("# phases considerees : "+phases)    
    
    print("# Repertoire des fichiers.txt : "+liste_chemin[ind])#fstat)
    print("# Fichier.data : "+jdd)
    print("# viscosite cinematique : "+str(nu))
    print("######################################")
    
    
    ## Preparation nom generique des fichiers utilises
    fstat = glob.glob(liste_chemin[ind]+"/diphasique_moyenne_spatiale_*.txt")
    if fstat == []:
        fstat = glob.glob(liste_chemin[ind]+"/monophasique_moyenne_spatiale_*.txt")
        print("Simulation monophasique")
    else :
        print("Simulation diphasique")


    # Tri des fichiers fstat par ordre croissant des seccondes. 
    # fstat.sort() --> Mauvais car place la seconde 10 juste apres la fin de la seconde 2, et pas 9
    fstat.sort(key=natural_keys)
    ################# Pour récuperer les numéros de colonne qu'il nous faut #########
    n_col_I = (2-1)
    n_col_uI = (3-1,4-1,5-1)
    n_col_uIv = (7-1,8-1,9-1)
    n_col_uuI = (11-1,12-1,13-1,14-1,15-1,16-1)
    n_col_dissip = (228-1)
    n_col_diujI = (72-1,73-1,74-1,75-1,76-1,77-1,78-1,79-1,80-1)
    #################################################################################
    t = []
    UL, VL, WL = [],[],[]
    UV, VV, WV = [],[],[]
    
    for f in fstat:
        t.append(float("."+f.split(".")[-2]))
        sec = float(f.split(".")[-3].split("_")[-1])
        t[-1]+=sec
        #
        if ((t[-1]>=t_deb) and (t[-1]<=t_fin)): 
            if (("liq" in phases) or ("all" in phases)):
                # Lecture
                I = np.loadtxt(f,usecols=(n_col_I)).T
                u_l = np.loadtxt(f,usecols=(n_col_uI[0])).T
                v_l = np.loadtxt(f,usecols=(n_col_uI[1])).T
                w_l = np.loadtxt(f,usecols=(n_col_uI[2])).T
                # Moyenne spatiale
                m_I = np.mean(I,axis=-1)
                m_u_l = np.mean(u_l,axis=-1)
                m_v_l = np.mean(v_l,axis=-1)
                m_w_l = np.mean(w_l,axis=-1)
                # Moyenne de phase
                m_u_ll = m_u_l/m_I  
                m_v_ll = m_v_l/m_I  
                m_w_ll = m_w_l/m_I  
                # Ajout à la liste
                UL.append(m_u_ll);VL.append(m_v_ll);WL.append(m_w_ll)
    
            if (("vap" in phases) or ("all" in phases)):
                # Lecture
                v = 1 - np.loadtxt(f,usecols=(n_col_I)).T
                u_v = np.loadtxt(f,usecols=(n_col_uIv[0])).T
                v_v = np.loadtxt(f,usecols=(n_col_uIv[1])).T
                w_v = np.loadtxt(f,usecols=(n_col_uIv[2])).T
                # Moyenne spatiale
                m_v = np.mean(v,axis=-1)
                m_u_v = np.mean(u_v,axis=-1)
                m_v_v = np.mean(v_v,axis=-1)
                m_w_v = np.mean(w_v,axis=-1)
                # Moyenne de phase
                m_u_vv = m_u_v/m_v
                m_v_vv = m_v_v/m_v
                m_w_vv = m_w_v/m_v
                # Ajout à la liste
                UV.append(m_u_vv);VV.append(m_v_vv);WV.append(m_w_vv)
        else:
            t.pop()
    
    
    
    ##################### Pour les  Plots #########################
    #Objectif : gérer l'affichage des ylabel en fonction des valeurs de UV
    if ("xyz_liq-vap" in composantes+"_"+phases):
        MUV,MUL,MVV,MVL,MWV,MWL = max(UV),max(UL),max(VV),max(VL),max(WV),max(WL)
        mUV,mUL,mVV,mVL,mWV,mWL = min(UV),min(UL),min(VV),min(VL),min(WV),min(WL)
        MU, mU = max([MUV,MUL,MVV,MVL,MWV,MWL]), min([mUV,mUL,mVV,mVL,mWV,mWL])
        
        coeff = np.floor(mt.log10(MU))
        ylabs = np.linspace(mU,MU,7)
        
        print("min(U),max(U),coeff,min(t),max(t)")
        print(MU,mU,coeff)
    # xlabs = np.floor(100*(np.linspace(min(t),max(t),7)))*10**(-2)
    
    ############################### Plots #########################
    print(t)
    dt = t[-1]-t[0]
    pas_pour_pente = int(len(t)/4)
    filtre_pour_smooth = np.ones(int(pas_pour_pente))/pas_pour_pente
    print("# t = "+str(t))
    print("# dt = "+str(dt))
    
    
    if (("x" in composantes) and ("liq" in phases)):
        d_ul = UL[-1]-UL[0]
        d_uldt = d_ul/dt; puiss = np.floor(mt.log10(abs(d_uldt)))
        print("# Derive u_x liquide : d_ul = "+str(d_ul))
        print("                       d_ul/dt ="+str(d_uldt))
        ax.plot(t, UL,'r',label=r"$\overline{U_x}^l$; $\partial_t \overline{U_x}^l = $"
                                +str(np.round(d_uldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                +" : "+str(liste_signes[ind]),
                           linestyle=linestyle2[ind],
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])
                           
        if ("temoin" in liste_signes[ind]) :                           
            ax_corr.plot(t, UL,'r',label=r"$\overline{U_x}^l$; $\partial_t \overline{U_x}^l = $"
                                        +str(np.round(d_uldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                        +" : "+str(liste_signes[ind]),
                                   linewidth=linewidth[ind],
                                   alpha=alpha[ind])
            ax_cp.plot(t, UL,'r',label=r"$\overline{U_x}^l$; $\partial_t \overline{U_x}^l = $"
                                      +str(np.round(d_uldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                      +" : "+str(liste_signes[ind]),
                                 linewidth=linewidth[ind],
                                 alpha=alpha[ind])
            UL_temoin_init = UL[0]
        
        # Pour une comparaison plus facile
        if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :        
            ax_corr.plot(t, UL-UL[0]+UL_temoin_init,'r',label=r"$\overline{U_x}^l$; $\partial_t \overline{U_x}^l = $"
                                                              +str(np.round(d_uldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                              +" : "+str(liste_signes[ind]+"*"),
                                                         linewidth=linewidth[ind],
                                                         alpha=alpha[ind])    
            ax_cp.plot(t, UL-UL[0]+UL_temoin_init,'r',label=r"$\overline{U_x}^l$; $\partial_t \overline{U_x}^l = $"
                                                           +str(np.round(d_uldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                           +" : "+str(liste_signes[ind]+"*"),
                                                      linewidth=linewidth[ind],
                                                      alpha=alpha[ind])    
        # Ajout des pentes
        pente_init, pente_fina = 0, 0
        # for i in range(int(len(t)/5)):
        #     pente_init += (UL[i]-UL[i+1]) / (t[i]-t[i+1])
        #     pente_fina += (UL[-i-1]-UL[-i]) / (t[-i-1]-t[-i])               
        pente_init = (-UL[0]+UL[pas_pour_pente]) / (-t[0]+t[pas_pour_pente])
        pente_fina = (UL[-1]-UL[-1-pas_pour_pente]) / (t[-1]-t[-1-pas_pour_pente])
        dt_init = -t[0]+t[pas_pour_pente]
        dt_fina = -t[-1-pas_pour_pente]+t[-1]
        ax_cp.plot([0,dt_init],[0.02,0.02+dt_init*pente_init],'r',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind],
                            label=r"$\partial_t \overline{U_x}^l_{init} = $"+str(np.round(pente_init*10**(-puiss),2))+r"$\times 10^{%s}$"%(str(puiss))+
                                  r"; $\partial_t \overline{U_x}^l_{final} = $"+str(np.round(pente_fina*10**(-puiss),2))+r"$\times 10^{%s}$"%(str(puiss)))
                            
        ax_cp.plot([5-dt_fina,5],[0.02-dt_fina*pente_fina,0.02],'r',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind])
        
        ax_pente.plot([0,dt_init],[0.02,0.02+dt_init*pente_init],'r',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind],
                            label=r"$\partial_t \overline{U_x}^l_{init} = $"+str(np.round(pente_init*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                  r"; $\partial_t \overline{U_x}^l_{final} = $"+str(np.round(pente_fina*10**(-puiss),2))+r"$\times 10^{%s}$"%(str(puiss)))
                            
        ax_pente.plot([5-dt_fina,5],[0.02-dt_fina*pente_fina,0.02],'r',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind])
        

        
    if (("y" in composantes) and ("liq" in phases)):
        print(linestyle_rev)
        d_vl = VL[-1]-VL[0]
        d_vldt = d_vl/dt; puiss = np.floor(mt.log10(abs(d_vldt)))
        print("# Derive u_y liquide : d_vl = "+str(d_vl))
        print("                       d_vl/dt ="+str(d_vldt))
        ax.plot(t, VL,'g',label=r"$\overline{U_y}^l$; $\partial_t \overline{U_y}^l = $"
                                +str(np.round(d_vldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           linestyle=linestyle2[ind],
                           alpha=alpha[ind])
                           
        # Pour une comparaison plus facile
        if ("temoin" in liste_signes[ind]) :                           
            ax_corr.plot(t, VL,'g',label=r"$\overline{U_y}^l$; $\partial_t \overline{U_y}^l = $"
                                        +str(np.round(d_vldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                        +" : "+str(liste_signes[ind]),
                                   linewidth=linewidth[ind],
                                   alpha=alpha[ind])
            VL_temoin_init = VL[0]
            ax_cp.plot(t, VL,'g',label=r"$\overline{U_y}^l$; $\partial_t \overline{U_y}^l = $"
                                       +str(np.round(d_vldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                       +" : "+str(liste_signes[ind]),
                                  linewidth=linewidth[ind],
                                  alpha=alpha[ind])
            
        if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :        
            ax_corr.plot(t, VL-VL[0]+VL_temoin_init,'g',label=r"$\overline{U_y}^l$; $\partial_t \overline{U_y}^l = $"
                                                             +str(np.round(d_vldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                             +" : "+str(liste_signes[ind]+"*"),
                                                        linewidth=linewidth[ind],
                                                        alpha=alpha[ind])
            ax_cp.plot(t, VL-VL[0]+VL_temoin_init,'g',label=r"$\overline{U_y}^l$; $\partial_t \overline{U_y}^l = $"
                                                           +str(np.round(d_vldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                           +" : "+str(liste_signes[ind]+"*"),
                                                      linewidth=linewidth[ind],
                                                      alpha=alpha[ind])
        # Ajout des pentes
        pente_init, pente_fina = 0, 0
        # for i in range(int(len(t)/5)):
        #     pente_init += (VL[i]-VL[i+1]) / (t[i]-t[i+1])
        #     pente_fina += (VL[-i-1]-VL[-i]) / (t[-i-1]-t[-i])               
        pente_init = (-VL[0]+VL[pas_pour_pente]) / (-t[0]+t[pas_pour_pente])
        pente_fina = (VL[-1]-VL[-1-pas_pour_pente]) / (t[-1]-t[-1-pas_pour_pente])               
        dt_init = -t[0]+t[pas_pour_pente]
        dt_fina = t[-1]-t[-1-pas_pour_pente] 
        ax_pente.plot([0,dt_init],[0.02,0.02+dt_init*pente_init],'g',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind],
                            label=r"$\partial_t \overline{U_x}^l_{init} = $"+str(np.round(pente_init*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                  r"; $\partial_t \overline{U_x}^l_{final} = $"+str(np.round(pente_fina*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss))
        ax_pente.plot([5-dt_fina,5],[0.02-dt_fina*pente_fina,0.02],'g',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind])
        ax_cp.plot([0,dt_init],[0.02,0.02+dt_init*pente_init],'g',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind],
                            label=r"$\partial_t \overline{U_x}^l_{init} = $"+str(np.round(pente_init*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                  r"; $\partial_t \overline{U_x}^l_{final} = $"+str(np.round(pente_fina*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss))
        ax_cp.plot([5-dt_fina,5],[0.02-dt_fina*pente_fina,0.02],'g',
                            linewidth=1,
                            linestyle="dashed",
                            alpha=alpha[ind])
                             
                                     
    if (("z" in composantes) and ("liq" in phases)):
        d_wl = WL[-1]-WL[0]
        d_wldt = d_wl/dt; puiss = np.floor(mt.log10(abs(d_wldt)))
        print("# Derive u_z liquide : d_wl = "+str(d_wl))
        print("                       d_wl/dt = "+str(d_wldt))    
        ax.plot(t, WL,'b',label=r"$\overline{U_z}^l$; $\partial_t \overline{U_z}^l = $"
                                +str(np.round(d_wldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                +" : "+str(liste_signes[ind]),
                          linewidth=linewidth[ind],
                          linestyle=linestyle2[ind],
                          alpha=alpha[ind])
                           
        # Pour une comparaison plus facile
        if ("temoin" in liste_signes[ind]) :                           
            ax_corr.plot(t, WL,'b',label=r"$\overline{U_z}^l$; $\partial_t \overline{U_z}^l = $"
                                         +str(np.round(d_wldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                         +" : "+str(liste_signes[ind]),
                                   linewidth=linewidth[ind],
                                   alpha=alpha[ind])   
            ax_cp.plot(t, WL,'b',label=r"$\overline{U_z}^l$; $\partial_t \overline{U_z}^l = $"
                                       +str(np.round(d_wldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                       +" : "+str(liste_signes[ind]),
                                 linewidth=linewidth[ind],
                                 alpha=alpha[ind])
            WL_temoin_init = WL[0]
        
        if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :        
            ax_corr.plot(t, WL-WL[0]+WL_temoin_init,'b',label=r"$\overline{U_z}^l$; $\partial_t \overline{U_z}^l = $"
                                                              +str(np.round(d_wldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                              +" : "+str(liste_signes[ind]+"*"),
                                                        linewidth=linewidth[ind],
                                                        alpha=alpha[ind])  
            ax_cp.plot(t, WL-WL[0]+WL_temoin_init,'b',label=r"$\overline{U_z}^l$; $\partial_t \overline{U_z}^l = $"
                                                            +str(np.round(d_wldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                            +" : "+str(liste_signes[ind]+"*"),
                                                      linewidth=linewidth[ind],
                                                      alpha=alpha[ind])                                                   
        # Ajout des pentes
        pente_init, pente_fina = 0, 0
        # for i in range(int(len(t)/5)):
        #     pente_init += (WL[i]-WL[i+1]) / (t[i]-t[i+1])
        #     pente_fina += (WL[-i-1]-WL[-i]) / (t[-i-1]-t[-i])               
        pente_init = (-WL[0]+WL[pas_pour_pente]) / (-t[0]+t[pas_pour_pente])
        pente_fina = (WL[-1]-WL[-1-pas_pour_pente]) / (t[-1]-t[-1-pas_pour_pente])               
        dt_init = -t[0]+t[pas_pour_pente] 
        dt_fina = t[-1]-t[-1-pas_pour_pente]
        ax_cp.plot([0,dt_init],[0.02,0.02+dt_init*pente_init],'b',
                         linewidth=1,
                         linestyle="dashed",
                         alpha=alpha[ind],
                         label=r"$\partial_t \overline{U_x}^l_{init} = $"+str(np.round(pente_init*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                  r"; $\partial_t \overline{U_x}^l_{final} = $"+str(np.round(pente_fina*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss))
                         
        ax_cp.plot([5-dt_fina,5],[0.02-dt_fina*pente_fina,0.02],'b',
                         linewidth=1,
                         linestyle="dashed",
                         alpha=alpha[ind])

        ax_pente.plot([0,dt_init],[0.02,0.02+dt_init*pente_init],'b',
                         linewidth=1,
                         linestyle="dashed",
                         alpha=alpha[ind],
                         label=r"$\partial_t \overline{U_x}^l_{init} = $"+str(np.round(pente_init*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                  r"; $\partial_t \overline{U_x}^l_{final} = $"+str(np.round(pente_fina*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss))
                         
        ax_pente.plot([5-dt_fina,5],[0.02-dt_fina*pente_fina,0.02],'b',
                         linewidth=1,
                         linestyle="dashed",
                         alpha=alpha[ind])
                           
    ############## VAPEUR
    ################### les figures avec corrections et pentes ne sont pas DU TOUT operationnelles    
    if (("x" in composantes) and ("vap" in phases)):
        ## RAW
        d_uv = UV[-1]-UV[0]
        d_uvdt = d_uv/dt; puiss = np.floor(mt.log10(abs(d_uvdt)))
        print("# Derive u_x vapeur : d_uv = "+str(d_uv))
        print("                       d_uv/dt = "+str(d_uvdt)) 
        ax.plot(t, UV,'r--',label=r"$\overline{U_x}^v$; $\partial_t \overline{U_x}^v = $"
                                  +str(np.round(d_uvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])
        ## SMOOTH
        d_uv_smt = np.convolve(UV,filtre_pour_smooth,mode='same')[-1]-np.convolve(UV,filtre_pour_smooth,mode='same')[0]
        d_uvdt_smt = d_uv_smt/dt; puiss_smt = np.floor(mt.log10(abs(d_uvdt_smt)))
        print("# smooth #")
        print("# Derive u_x vapeur : d_uv = "+str(d_uv_smt))
        print("                       d_uv/dt = "+str(d_uvdt_smt)) 
        ax_smt.plot(t, np.convolve(UV,filtre_pour_smooth,mode='same'),'r--',label=r"$\widetilde{\overline{U_x}^v}$; $\partial_t \widetilde{\overline{U_x}^v} = $"
                                  +str(np.round(d_uvdt_smt*10**(-puiss_smt),2))+r"$\times 10^{%s}$"%(puiss_smt)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])

    if ("temoin" in liste_signes[ind]) :    
        UV_temoin_init = UV[0]
        
      #  if ("temoin" in liste_signes[ind]) :                           
      #      ax_corr.plot(t, UV,'r--',label=r"$\overline{U_x}^v$; $\partial_t \overline{U_x}^v = $"
      #                              +str(np.round(d_uvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
      #                              +" : "+str(liste_signes[ind]),
      #                      linewidth=linewidth[ind],
      #                      alpha=alpha[ind])

      #  if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :        
      #      ax_corr.plot(t, UV-UV[0]+UV_temoin_init,'r',label=r"$\overline{U_x}^l$; $\partial_t \overline{U_x}^l = $"
      #                                                +str(np.round(d_uldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
      #                                                +" : "+str(liste_signes[ind]+"*"),
      #                                             linewidth=linewidth[ind],
      #                                             alpha=alpha[ind])  
                                                   
        # Ajout des pentes (en l'etat du 07.07.21, si on ne donne pas temoin dans la liste des signes, la figure de correction n'aura que les panetes initiales et finales) NON OPERATIONEL POUR LA PHASE VAPEUR
      #  pente_init, pente_fina = 0, 0
      #  for i in range(int(len(t)/5)):
      #      pente_init += (UL[i]-UL[i+1]) / (t[i]-t[i+1])
      #      pente_fina += (UL[-i-1]-UL[-i]) / (t[-i-1]-t[-i])               
      #  dt_init = -t[0]+t[int(len(t)/4)]
      #  dt_fina = -t[-int(len(t)/4)]+t[-1]
      #  ax_corr.plot([t[0],t[int(len(t)/4)]],[UL[0],UL[0]+dt_init*pente_init],'r',
      #                   linewidth=1,
      #                   linestyle="dashed",
      #                   alpha=alpha[ind],
      #                   label="pente init : "+str(pente_init)+"; pente finale : "+str(pente_fina))
      #                   
      #  ax_corr.plot([t[-int(len(t)/4)],t[-1]],[UL[-1]-dt_fina*pente_fina,UL[-1]],'r',
      #                   linewidth=1,
      #                   linestyle="dashed",
      #                   alpha=alpha[ind])
        
    if (("y" in composantes) and ("vap" in phases)):
        ## RAW
        d_vv = VV[-1]-VV[0]
        d_vvdt = d_vv/dt; puiss = np.floor(mt.log10(abs(d_vvdt)))
        print("# Derive u_y vapeur : d_vV = "+str(d_vv))
        print("                      d_vV/dt = "+str(d_vvdt))    
        ax.plot(t, VV,'g--',label=r"$\overline{U_y}^v$; $\partial_t \overline{U_y}^v = $"
                                  +str(np.round(d_vvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])

        if ("temoin" in liste_signes[ind]) :                           
            ax.plot(t, VV,'g--',label=r"$\overline{U_y}^v$; $\partial_t \overline{U_y}^v = $"
                                      +str(np.round(d_vvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                    +" : "+str(liste_signes[ind]),
                               linewidth=linewidth[ind],
                               alpha=alpha[ind])
            VL_temoin_init = VL[0]

        ## SMOOTH
        d_vv_smt = np.convolve(VV,filtre_pour_smooth,mode='same')[-1]-np.convolve(VV,filtre_pour_smooth,mode='same')[0]
        d_vvdt_smt = d_vv_smt/dt; puiss_smt = np.floor(mt.log10(abs(d_vvdt_smt)))
        print("# smooth #")
        print("# Derive u_y vapeur : d_wl = "+str(d_vv_smt))
        print("                       d_wl/dt = "+str(d_vvdt_smt)) 
        ax_smt.plot(t, np.convolve(VV,filtre_pour_smooth,mode='same'),'g--',label=r"$\widetilde{\overline{U_y}^v}$; $\partial_t \widetilde{\overline{U_y}^v} = $"
                                  +str(np.round(d_vvdt_smt*10**(-puiss_smt),2))+r"$\times 10^{%s}$"%(puiss_smt)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])
        
        if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :        
            ax_corr.plot(t, VL-VL[0]+VL_temoin_init,'g',label=r"$\overline{U_y}^l$; $\partial_t \overline{U_y}^l = $"
                                                       +str(np.round(d_vldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                                       +" : "+str(liste_signes[ind]+"*"),
                                                  linewidth=linewidth[ind],
                                                  alpha=alpha[ind])
                                                  
        # Ajout des pentes (en l'etat du 07.07.21, si on ne donne pas temoin dans la liste des signes, la figure de correction n'aura que les panetes initiales et finales)
      #  pente_init, pente_fina = 0, 0
      #  for i in range(int(len(t)/5)):
      #      pente_init += (VL[i]-VL[i+1]) / (t[i]-t[i+1])
      #      pente_fina += (VL[-i-1]-VL[-i]) / (t[-i-1]-t[-i])               
      #  dt_init = -t[0]+t[int(len(t)/4)]
      #  dt_fina = -t[-int(len(t)/4)]+t[-1]
      #  ax_corr.plot([t[0],t[int(len(t)/4)]],[VL[0],VL[0]+dt_init*pente_init],'g',
      #                      linewidth=1,
      #                      linestyle="dashed",
      #                      alpha=alpha[ind],
      #                      label="pente init : "+str(pente_init)+"; pente finale : "+str(pente_fina))
      #                      
      #  ax_corr.plot([t[-int(len(t)/4)],t[-1]],[VL[-1]-dt_fina*pente_fina,VL[-1]],'g',
      #                      linewidth=1,
      #                      linestyle="dashed",
      #                      alpha=alpha[ind])
        
    if (("z" in composantes) and ("vap" in phases)):
        ## RAW
        d_wv = WV[-1]-WV[0]
        d_wvdt = d_wv/dt; puiss = np.floor(mt.log10(abs(d_wvdt)))
        print("# Derive u_z vapeur : d_wv/dt = "+str(d_wv))
        print("                       d_wv/dt = "+str(d_wvdt)) 
        ax.plot(t, WV,'b--',label=r"$\overline{U_z}^v$; $\partial_t \overline{U_z}^v = $"
                                  +str(np.round(d_wvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])

        ## SMOOTH
        d_wv_smt = np.convolve(WV,filtre_pour_smooth,mode='same')[-1]-np.convolve(WV,filtre_pour_smooth,mode='same')[0]
        d_wvdt_smt = d_wv_smt/dt; puiss_smt = np.floor(mt.log10(abs(d_wvdt_smt)))
        print("# smooth #")
        print("# Derive u_z vapeur : d_wv = "+str(d_wv_smt))
        print("                       d_wv/dt = "+str(d_wvdt_smt)) 
        ax_smt.plot(t, np.convolve(WV,filtre_pour_smooth,mode='same'),'b--',label=r"$\widetilde{\overline{U_z}^v}$; $\partial_t \widetilde{\overline{U_z}^v} = $"
                                  +str(np.round(d_wvdt_smt*10**(-puiss_smt),2))+r"$\times 10^{%s}$"%(puiss_smt)
                                +" : "+str(liste_signes[ind]),
                           linewidth=linewidth[ind],
                           alpha=alpha[ind])

        if ("temoin" in liste_signes[ind]) :                           
            ax.plot(t, WV,'b--',label=r"$\overline{U_z}^v$; $\partial_t \overline{U_z}^v = $"
                                      +str(np.round(d_wvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
                                    +" : "+str(liste_signes[ind]),
                               linewidth=linewidth[ind],
                               alpha=alpha[ind])
            WV_temoin_init = WV[0]

      #  if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :        
      #      ax_corr.plot(t, WL-WL[0]+WL_temoin_init,'b',label=r"$\overline{U_z}^l$; $\partial_t \overline{U_z}^l = $"
      #                                                  +str(np.round(d_wldt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)
      #                                                +" : "+str(liste_signes[ind]+"*"),
      #                                           linewidth=linewidth[ind],
      #                                           alpha=alpha[ind])
      #                                           
      #      # Ajout des pentes
      #      pente_init, pente_fina = 0, 0
      #      for i in range(int(len(t)/5)):
      #          pente_init += (WL[i]-WL[i+1]) / (t[i]-t[i+1])
      #          pente_fina += (WL[-i-1]-WL[-i]) / (t[-i-1]-t[-i])               
      #      dt_init = -t[0]+t[int(len(t)/4)]
      #      dt_fina = -t[-int(len(t)/4)]+t[-1]
      #      ax_corr.plot([t[0],t[int(len(t)/4)]],[WL[0],WL[0]+dt_init*pente_init],'b',
      #                       linewidth=1,
      #                       linestyle="dashed",
      #                       alpha=alpha[ind])
      #      ax_corr.plot([t[-int(len(t)/4)],t[-1]],[WL[-1]-dt_fina*pente_fina,WL[-1]],'b',
      #                       linewidth=1,
      #                       linestyle="dashed",
      #                       alpha=alpha[ind])
                             
        
    # ax.text(0.99*(t[-1]+t[0])/2, 0.8*UL, r'$\langle{K_F} \rangle = $ '+str(np.round(mUV*10**(-mUV_coeff),5))+r'$\cdot 10^{%s}$'%(mUV_coeff), fontsize=18)
    ax.set_xlabel('temps (s)',fontsize=18)
    ax_smt.set_xlabel('temps (s)',fontsize=18)
    ax_corr.set_xlabel('temps (s)',fontsize=18)
    ax_cp.set_xlabel('temps (s)',fontsize=18)
    ax_pente.set_xlabel('temps (s)',fontsize=18)
    # ax.set_xticks(xlabs)
    # ax.set_xticklabels(np.round(xlabs,2),fontsize=18)
    ax.axhline(y=0, color='k',linewidth=0.5)
    ax_smt.axhline(y=0, color='k',linewidth=0.5)
    ax_corr.axhline(y=0, color='k',linewidth=0.5)
    ax_cp.axhline(y=0, color='k',linewidth=0.5)
    # ax_pente.axhline(y=0, color='k',linewidth=0.5)
    # if ("xyz_liq-vap" in composantes+"_"+phases):
        # ax.set_ylabel(r'$\overline{U_i}^k \times 10^{%s}$'%(str(coeff)),fontsize=18)
        # ax.set_yticks(ylabs)
        # ax.set_yticklabels(np.round(ylabs*10**(-coeff),3),fontsize=18)
    # else:
    ax.set_ylabel(r'$\overline{U_i}^k (m.s^{-1})$',fontsize=18)
    ax_smt.set_ylabel(r'$\overline{U_i}^k (m.s^{-1})$',fontsize=18)
    ax_corr.set_ylabel(r'$\overline{U_i}^k (m.s^{-1})$ corrigees',fontsize=18)
    ax_cp.set_ylabel(r'$\overline{U_i}^k (m.s^{-1})$ corrigees',fontsize=18)
    ax_pente.set_ylabel(r'$\overline{U_i}^k (m.s^{-1})$ corrigees',fontsize=18)
    

ax.tick_params(labelsize=18)
ax_smt.tick_params(labelsize=18)
ax_corr.tick_params(labelsize=18)
ax_cp.tick_params(labelsize=18)
ax_pente.tick_params(labelsize=18)
ax.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5))
ax_smt.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5))
ax_corr.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5))
ax_cp.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5))
ax_pente.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5))
fig.tight_layout()
fig_smt.tight_layout()
fig_corr.tight_layout()
fig_cp.tight_layout()
fig_pente.tight_layout()
fig.savefig(figname+"_"+composantes+"_"+phases+".png")
fig_smt.savefig(figname+"_"+composantes+"_"+phases+"_smt"+".png")
fig_corr.savefig(figname+"_"+composantes+"_"+phases+"_corr"+".png")
fig_cp.savefig(figname+"_"+composantes+"_"+phases+"_cp"+".png")
fig_pente.savefig(figname+"_"+composantes+"_"+phases+"_pente"+".png")
fig.savefig(figname+"_"+composantes+"_"+phases+".pdf")
fig_corr.savefig(figname+"_"+composantes+"_"+phases+"_corr"+".pdf")
fig_cp.savefig(figname+"_"+composantes+"_"+phases+"_cp"+".pdf")
fig_pente.savefig(figname+"_"+composantes+"_"+phases+"_pente"+".pdf")
