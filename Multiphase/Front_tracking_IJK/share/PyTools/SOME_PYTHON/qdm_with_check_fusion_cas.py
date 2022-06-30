# -*- coding: utf8

import os
import numpy
import subprocess
import scipy
from scipy import signal
# import pandas
import glob
import sys


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
# matplotlib.use("Agg")
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

import math as mt

"""
    Fonctionne que pour check_etapes_et_termes.out .
"""

########################################################################
# mode d'emploi : python qdm_with_check_fusion_cas.py nom_fig liste_chemin_out liste_signes   directions  N_smooth [temps0,temps1] correction_listes
#                            0             1           2               3             4         5            6
# /!\ 3 : xyz ou xy ou x ou ... . all revient a xyz
# /!\ 4 : epaisseur du filtre pour adoucir les coubres 
#         (100 est ok pour le jdd de spectral_UO_enable_diphasique)
# /!\ 5 : [temps0,temps1] -> releve les donnees des outs entre les 2 temps
# /!\ 5 : [temps0,temps0] -> releve la donnee des outs a l'instant le plus proche de temps0
#         [-1,-1] -> releve tous les temps disponibles
# Exemple sur gr262753@is243146 :
# /volatile/RANDOM_TRIPERIO/SOME_PYTHON$ python SOME_PYTHON/quantite_mvt_old_bilan_fusion_cas.py fig_qdm_mix [petit_a015/COARSE/RUN04/OUT/,petit_a030/COARSE/RUN04/OUT/,petit_a060/COARSE/RUN06/OUT/,petit_a120/COARSE/RUN06/OUT/] [1.5,3.0,6.0,12] yz 10 [-1,-1]
########################################################################

message = "#############################################\n\
Tracons la QDM et les forces interfaciales au cours du temps \n\
##############################################"
print(message)

## Lecture des donnees utilisateur
nom_figure = sys.argv[1]
liste_signes = sys.argv[3].rstrip("]").lstrip("[").split(",")
directions = sys.argv[4]
N_smooth = float(sys.argv[5])
if directions == "all" : directions = "xyz"
########################################################################

# Travail sur les directions a considerer ##############################
liste_directions = []
if "x" in directions.lower() :
    liste_directions.append("x")
if "y" in directions.lower() :
    liste_directions.append("y")
if "z" in directions.lower() :
    liste_directions.append("z")


## Creation des chemins, recuperation du temps
liste_chemin_out = sys.argv[2].rstrip("]").lstrip("[").split(",")
good_liste_out = []
BIG_dix = int(sys.argv[7])
print("BIG_dix // correction : ", BIG_dix)
if BIG_dix:
  for out in liste_chemin_out:
    if (("N01" in out) or ("SIMU_120" in out)) :
        good_liste_out.append(out)
  liste_chemin_out = good_liste_out

### A SUPPRIMER AUSSITOT QUE LA SIMULAITON CONCERNEE A FINI DE TOURNER
good_liste_out = []
for out in liste_chemin_out:
  print("out : ",out)
  out = out.replace("00_1PROCY/RUN06/OU","00_1PROCY/RUN04/OU")  
  print("out : ",out)
  print("# mmm")
  good_liste_out.append(out)
liste_chemin_out = good_liste_out
# BIG_SIMU/RE_ISOLEE_143_G00_1PROCY/RUN03
###
linewidth = [1.,2.,4.,4.]
alpha = [1.0,0.6,0.4,0.25]

########################################################################
# Pierre de Rosette pour savoir quelle colonne correspond à quelle grandeur
#              car dans le cas d'une reprise il se peut que l'on perde cette information
#              meme dans le cas d'un calcul initial, il a de fortes chances pour que ce
#              ne soit pas "bien" fait.
########################################################################
rosette = { "it":0,
        "t":1,
        "rk":2,
        "ru_av_pred"         : {"x":3,  "y":4,  "z":5},   # rho_u_euler_av_pred
        "rdu_ap_pred"        : {"x":6,  "y":7,  "z":8},   # rho_du_euler_ap_pred
        "ru_ap_proj"         : {"x":9,  "y":10, "z":11},  # rho_u_euler_av_prediction
        "rdu_ap_proj"        : {"x":12, "y":13, "z":14},  # rho_du_euler_ap_prediction
        "ru_av_rmi"          : {"x":15, "y":16, "z":17},  # rho_u_euler_av_rho_mu_ind
        "ru_ap_rmi"          : {"x":18, "y":19, "z":20},  # rho_u_euler_ap_rho_mu_ind
        "u_ap_rmi"           : {"x":21, "y":22, "z":23},  # u_euler_ap_rho_mu_ind
        "t_intf"             : {"x":24, "y":25, "z":26},  # terme_interfaces               /!\/!\ il est nul celui-là, on le remplit que sous certaines conditions... ne le regardons pas
        "t_conv"             : {"x":27, "y":28, "z":29},  # terme_convection : defaut -> Moyenne_spatiale{ \rho S [u_i u_j] }   @\rho du/dt
        "t_diff"             : {"x":30, "y":31, "z":32},  # terme_diffusion : defaut -> Moyenne_spatiale{ S [\mu s_{ij}] } *    @\rho du/dt
        "t_pr_2"             : {"x":33, "y":34, "z":35},  # terme_pression_bis : Moyenne_spatiale{ grad(p*) }                   @\rho du/dt
        "t_pr_3"             : {"x":36, "y":37, "z":38},  # terme_pression_bis : Moyenne_spatiale{ 1/rho * grad(p*) }           @     du/dt

        "t_intf_bf_ms_bis"   : {"x":39, "y":40, "z":41},  # terme_interfaces_bf_mass_solveur_bis : Moyenne_spatiale{ dv-tc-td }                 @\rho du/dt
        "t_intf_bf_ms"       : {"x":42, "y":43, "z":44},  # terme_interfaces_bf_mass_solver : Moyenne_spatiale{ terme_source_interfaces_ns_ }   @\rho du/dt
        "t_intf_af_ms"       : {"x":45, "y":46, "z":47},  # terme_interfaces_af_mass_solver : Moyenne_spatiale{ terme_source_interfaces_ns_ }   @      du/dt (le mass_solver a ete applique)
        "pres_ap_proj"       : {"x":48, "y":48, "z":48},  # pression_ap_proj, ce n'est qu'UN double
        "t_conv_mass_sol"    : {"x":52, "y":53, "z":54},  # terme_moyen_convection_mass_solver_ : Moyenne_spatiale{ 1/rho * C }  @      du/dt
        "t_diff_mass_sol"    : {"x":55, "y":56, "z":57}}  # terme_moyen_diffusion_mass_solver_  : Moyenne_spatiale{ 1/rho * D }  @      du/dt

########################################################################
########################################################################

color = {"x":"r","y":"g","z":"b"}
linestyle = ["solid", "dashed", "dotted", "dashdot"]
print(" len(liste_chemin_out) : ",len(liste_chemin_out))
if len(liste_chemin_out)>1:
    linestyle_rev = linestyle 
    linestyle_rev.reverse() 
# linestyle2 : loosely dotted; loosely dashed; loosely dashdotted; solid (ofc
linestyle2 = [(0, (1,10)), (0, (5, 10)), (0, (3,10,1,10)), "solid"]
if len(liste_chemin_out)==1:
    linestyle2.reverse()
########################################################################

# Valeurs de references, propres au probleme
db = 1e-03
S = 4*numpy.pi*(db/2.)**2

## Sorties #############################################
### unsmoothed
fig_rv, (ax_rv) = plt.subplots(1,1,figsize=(14,7))
fig_rv_corr, (ax_rv_corr) = plt.subplots(1,1,figsize=(14,7))
fig_rv_pente, (ax_rv_pente) = plt.subplots(1,1,figsize=(14,7))
fig_rv_cp, (ax_rv_cp) = plt.subplots(1,1,figsize=(14,7))

fig_rv.suptitle(r"$\rho u$",fontsize=24)
fig_rv_corr.suptitle(r"$\rho u$, corrigees",fontsize=24)
fig_rv_pente.suptitle(r"$\rho u$, pentes",fontsize=24)
fig_rv_cp.suptitle(r"$\rho u$, corrigees et pentes",fontsize=24)

### smoothed
filtre = numpy.ones(int(N_smooth))/N_smooth
fig_smt_rv,    (ax_smt_rv) = plt.subplots(1,1,figsize=(14,7))
fig_smt_rv_corr,    (ax_smt_rv_corr) = plt.subplots(1,1,figsize=(14,7))
fig_smt_rv_pente,    (ax_smt_rv_pente) = plt.subplots(1,1,figsize=(14,7))
fig_smt_rv_cp,    (ax_smt_rv_cp) = plt.subplots(1,1,figsize=(14,7))

fig_smt_rv.suptitle(r"$\rho u$",fontsize=24)

rv_temoin_init = [0,0,0]
for ind, chemin_out in enumerate(liste_chemin_out):
    #chemin_out = chemin_out+"*bilan_qdm.out"
    chemin_out = chemin_out+"*check_etapes_et_termes.out"
    print("fichier out : "+chemin_out)
    t_deb, t_fin = float(sys.argv[6].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[6].rstrip("]").lstrip("[").split(",")[-1])
    print("# t_deb = %s \n# t_fin = %s"%(str(t_deb),str(t_fin)))
    
    print(numpy.loadtxt(glob.glob(chemin_out)[-1]).shape)
    temps = numpy.loadtxt(glob.glob(chemin_out)[-1])[:,1]
    ########################################################################
    
    
    ## Lecture, selection, rangement et traitement des donnees
    # chemin_out = "/volatile/RANDOM_TRIPERIO/a030/RE_ISOLE_143/RUN05/OUT/"
    
    bilan_qdm = numpy.loadtxt(glob.glob(chemin_out)[-1])#[rosette["ru_ap_rmi"]["x"]][it_deb:it_fin,1:]
    # print(bilan_qdm[rosette["tau_wa"]["x"]])
    
    ########################################################################
    
    # Travail sur les iterations ###########################################
    ## Recuperation des iterations a considerer dans les fichiers out
    # Si es instants renseignés ne sont pas inclus dans la plage de releves
    #   -> l'instant initial est le premier instant de releve. L'instant final est le dernier instant de releve
    if t_deb<=temps[0] : t_deb = temps[0]; print("Debut des releves : "+str(t_deb))
    if (t_fin<=0 or t_fin>temps[-1]) : t_fin = temps[-1]; print("Fin des releves : "+str(t_fin))
    # On associe les instants donnes a l'iteration correspondante
    it_deb = numpy.argwhere(numpy.abs(temps-t_deb)==numpy.min(numpy.abs(temps-t_deb)))
    it_fin = numpy.argwhere(numpy.abs(temps-t_fin)==numpy.min(numpy.abs(temps-t_fin)))
    # Si les instants donnes sont exactement a mi-temsp entre deux instants de releve :
    #    -> on choisi la plus petite valeur de debut et la plus grande valeur de fin
    if len(it_deb)>1:
        it_deb = it_deb[0]
    if len(it_fin)>1:
        it_fin = it_fin[-1]
    # Si on a donne le meme instant de debut et de fin :
    #   -> c'est comme si on avait donne un interval de longeur 1
    # if it_fin==it_deb : it_fin +=1
    
    it_deb, it_fin = int(it_deb), int(it_fin)
    print(it_deb, it_fin)
    temps = temps[it_deb:it_fin+1]
    DT = temps[-1] - temps[0]
        
    for i, direction in enumerate(liste_directions):

        ### NOT SMOOTHENED PLOTS ###########################################
                         
        # sum{alpha * rho * u} = rho u : qdm
        rv = bilan_qdm[it_deb:it_fin+1,[rosette["ru_ap_rmi"][direction]]]
        d_rv = rv[-1] - rv[0]
        d_rvdt = (d_rv/DT)[0]; 
        if d_rvdt==0:
            puiss=0
        else:
            puiss = numpy.floor(mt.log10(abs(d_rvdt)))
        ax_rv.plot(temps,rv, color=color[direction],
                   linewidth=linewidth[ind],
                   linestyle=linestyle2[ind],
                   alpha=alpha[ind], 
                   label=direction+
                         r" $\overline{ \rho u}; \partial_t (\overline{ \rho u}) = $"+
                         str(numpy.round(d_rvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                         " : "+liste_signes[ind])
        print(liste_signes[ind]+" : "+str(rv[0]))
        # Pour une comparaison plus facile : 
        if (not("temoin" in liste_signes[ind]) and ("temoin" in liste_signes)) :
             ## superposition du point initial avec celui du temoin     
             rv_corr = rv-rv[0]+rv_temoin_init[i]
             print("rv_temoin_init : ",rv_temoin_init)
             pente_init, pente_finale = 0, 0
             ax_rv_corr.plot(temps,rv_corr, color=color[direction],
                             linewidth=linewidth[ind],
                             linestyle="solid",
                             alpha=alpha[ind],
                             label=direction+
                                   r" $\overline{ \rho u}; \partial_t (\overline{ \rho u}) = $"+
                                   str(numpy.round(d_rvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                   " : "+liste_signes[ind]+
                                   "*")
             ax_rv_cp.plot(temps,rv_corr, color=color[direction],
                             linewidth=linewidth[ind],
                             linestyle="solid",
                             alpha=alpha[ind],
                             label=direction+
                                   r" $\overline{ \rho u}; \partial_t (\overline{ \rho u}) = $"+
                                   str(numpy.round(d_rvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                   " : "+liste_signes[ind]+
                                   "*")
             ## Pentes initiale et finale
             pas_pour_pente = int(len(temps)/4)
             # for i in range(int(len(temps)/5) :
             #     pente_init += (rv[i]-rv[i+1]) / (temps[i]-temps[i+1])
             #     pente_finale += (rv[-i-1]-rv[-i]) / (temps[-i-1]-temps[-i])
             pente_init = (-rv[0]+rv[pas_pour_pente]) / (-temps[0]+temps[pas_pour_pente])
             pente_fina = (rv[-1]-rv[-pas_pour_pente]) / (temps[-1]-temps[-pas_pour_pente])
             dt_init = -temps[0]+temps[pas_pour_pente]
             dt_fina = -temps[-pas_pour_pente]+temps[-1]
             ax_rv_pente.plot([0,dt_init],[5,5+dt_init*pente_init], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind],
                             label=r"$\partial_t (\overline{ \rho u})_{init} = $"+str(numpy.round(pente_init[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+"; "+
                                   r"$\partial_t (\overline{ \rho u})_{fina} = $"+str(numpy.round(pente_fina[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)) 
             ax_rv_pente.plot([5-dt_fina,5],[5-dt_fina*pente_fina,5], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind])#,
                             #label="pente finale : "+str(pente_fina[0]))
             ax_rv_cp.plot([0,dt_init],[5,5+dt_init*pente_init], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind],
                             label=r"$\partial_t (\overline{ \rho u})_{init} = $"+str(numpy.round(pente_init[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+"; "+
                                   r"$\partial_t (\overline{ \rho u})_{fina} = $"+str(numpy.round(pente_fina[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)) 
             ax_rv_cp.plot([5-dt_fina,5],[5-dt_fina*pente_fina,5], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind])#,
                             #label="pente finale : "+str(pente_fina[0]))
                         
                         
        # Pour une comparaison plus facile, temoin
        if ("temoin" in liste_signes[ind]):
             rv_temoin_init[i] = rv[0]
             pente_init_temoin, pente_finale_temoin = 0, 0
             ax_rv_corr.plot(temps,rv, color=color[direction],
                             linewidth=linewidth[ind],
                             linestyle="solid",
                             alpha=alpha[ind],
                             label=direction+
                                   r" $\overline{ \rho u}; \partial_t (\overline{ \rho u}) = $"+
                                   str(numpy.round(d_rvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                   " : "+liste_signes[ind]+
                                   "*")
             ax_rv_cp.plot(temps,rv, color=color[direction],
                             linewidth=linewidth[ind],
                             linestyle="solid",
                             alpha=alpha[ind],
                             label=direction+
                                   r" $\overline{ \rho u}; \partial_t (\overline{ \rho u}) = $"+
                                   str(numpy.round(d_rvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                                   " : "+liste_signes[ind]+
                                   "*")
             ## Pentes initiale et finale
             pas_pour_pente = int(len(temps)/4)
             # for i in range(int(len(temps)/5) :
             #     pente_init_temoin += (rv[i]-rv[i+1]) / (temps[i]-temps[i+1])
             #     pente_init_temoin += (rv[i]-rv[i+1]) / (temps[i]-temps[i+1])
             pente_init_temoin = (-rv[0]+rv[pas_pour_pente]) / (-temps[0]+temps[pas_pour_pente]) 
             pente_fina_temoin = (rv[-1]-rv[-pas_pour_pente]) / (temps[-1]-temps[-pas_pour_pente]) 
             dt_init = -temps[0]+temps[pas_pour_pente]
             dt_fina = -temps[-pas_pour_pente]+temps[-1]
             ax_rv_pente.plot([0,dt_init],[5,5+dt_init*pente_init_temoin], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind],
                             label=r"$\partial_t (\overline{ \rho u})_{init} = $ : "+str(numpy.round(pente_init_temoin[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+"; "+
                                   r"$\partial_t (\overline{ \rho u})_{fina} = $ : "+str(numpy.round(pente_fina_temoin[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss))
 
             ax_rv_pente.plot([5-dt_fina,5],[5-dt_fina*pente_finale_temoin,5], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind])#,
                             #label="pente finale : "+str(pente_init_temoin))             
             ax_rv_cp.plot([0,dt_init],[5,5+dt_init*pente_init_temoin], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind],
                             label=r"$\partial_t (\overline{ \rho u})_{init} = $ : "+str(numpy.round(pente_init_temoin[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+"; "+
                                   r"$\partial_t (\overline{ \rho u})_{fina} = $ : "+str(numpy.round(pente_fina_temoin[0]*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss))
 
             ax_rv_cp.plot([5-dt_fina,5],[5-dt_fina*pente_finale_temoin,5], color=color[direction],
                             linewidth=1,
                             linestyle="dashed",
                             alpha=alpha[ind])#,
                             #label="pente finale : "+str(pente_init_temoin))             
      

        # pour les autres termes --> voir : check_etapes_et_termes
        
        # Bilan
        # --> voir check_etapes_et_termes

        ### FIN ### NOT SMOOTHENED PLOTS ###############################
    
    
    
        #### SMOOTHENED PLOTS #######################################
        # Force interfaciales
        # rho v
        rv = numpy.convolve(bilan_qdm[it_deb:it_fin+1,[rosette["ru_ap_rmi"][direction]]][:,0],filtre,mode='same')
        d_rv = rv[-1] - rv[0]
        d_rvdt = d_rv/DT;
        if d_rvdt==0:
            puiss=0
        else:
            puiss = numpy.floor(mt.log10(abs(d_rvdt)))
        ax_smt_rv.plot(temps,rv, color=color[direction],
                       linewidth=linewidth[ind],
                       alpha=alpha[ind], 
                       label=direction+
                             r" $\overline{ \rho u}; \partial_t (\overline{ \rho u}) = $"+
                             str(numpy.round(d_rvdt*10**(-puiss),2))+r"$\times 10^{%s}$"%(puiss)+
                             " : "+str(liste_signes[ind]))
    
        # Convection, Diffusion et pression
        # --> voir check_etapes_et_termes
        
        # Bilan
        # --> voir check_etapes_et_termes
        #### FIN ## SMOOTHENED PLOTS #######################################   
    
    ### UNSMOOTH #########################################################
    ax_rv.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5));ax_rv.grid(True,'minor')
    ax_rv_corr.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5));ax_rv_corr.grid(True,'minor')
    ax_rv_pente.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5));ax_rv_pente.grid(True,'minor')
    ax_rv_cp.legend(fontsize=18,loc='center left', bbox_to_anchor=(1.04, 0.5));ax_rv_cp.grid(True,'minor')
    plt.tight_layout()

    ax_rv.set_xlabel("temps (s)",fontsize=24)
    ax_rv_corr.set_xlabel("temps (s)",fontsize=24)
    ax_rv_pente.set_xlabel("temps (s)",fontsize=24)
    ax_rv_cp.set_xlabel("temps (s)",fontsize=24)
    
    ax_rv.tick_params(labelsize=24)
    ax_rv_corr.tick_params(labelsize=24)
    ax_rv_pente.tick_params(labelsize=24)
    ax_rv_cp.tick_params(labelsize=24)
    
    ax_rv.set_ylabel("Quantite de mouvement",fontsize=24)
    ax_rv_corr.set_ylabel("Quantite de mouvement",fontsize=24)
    ax_rv_pente.set_ylabel("Quantite de mouvement",fontsize=24)
    ax_rv_cp.set_ylabel("Quantite de mouvement",fontsize=24)
    ### FIN # UNSMOOTH ###################################################

    ### SMOOTHED #######################################################
    ax_smt_rv.legend(fontsize=18);ax_smt_rv.grid(True,'minor')
    
    ax_smt_rv.set_xlabel("temps (s)",fontsize=24)
    
    ax_smt_rv.tick_params(labelsize=24)
    
    ax_smt_rv.set_ylabel("Quantite de mouvement",fontsize=24)
    ### FIN # SMOOTHENED ###############################################

ax_rv.axhline(y=0, color='k',linewidth=0.5)
ax_rv_corr.axhline(y=0, color='k',linewidth=0.5)
ax_rv_pente.axhline(y=0, color='k',linewidth=0.5)
ax_rv_cp.axhline(y=0, color='k',linewidth=0.5)
### UNSMOOTHED #########################################################
fig_rv.legend(loc='center left', bbox_to_anchor=(1.04, 0.5));fig_rv.tight_layout()
fig_rv_corr.legend(loc='center left', bbox_to_anchor=(1.04, 0.5));fig_rv_corr.tight_layout()
fig_rv_pente.legend(loc='center left', bbox_to_anchor=(1.04, 0.5));fig_rv_pente.tight_layout()
fig_rv_cp.legend(loc='center left', bbox_to_anchor=(1.04, 0.5));fig_rv_cp.tight_layout()

# plt.show()

fig_rv.savefig(nom_figure+"_rv_"+directions+".png")
fig_rv_corr.savefig(nom_figure+"_rv_corr_"+directions+".png")
fig_rv_pente.savefig(nom_figure+"_rv_pente_"+directions+".png")
fig_rv_cp.savefig(nom_figure+"_rv_cp_"+directions+".png")
### SMOOTHED ###########################################################
fig_smt_rv.savefig(nom_figure+"_rv_smooth_"+directions+".png")
