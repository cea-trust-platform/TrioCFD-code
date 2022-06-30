# -*- coding: utf8

import os
import numpy
import subprocess
import scipy
from scipy import signal
import pandas
import glob
import sys


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
# matplotlib.use("Agg")
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

"""
    Ne fonctionne pas encore. Travailler sur la fiche qui se trouve ici :
    /volatile/EQUATION_ENERGIE/IJK_EXTENSION_3_DIPHASIQUE_CDP
"""

########################################################################
# mode d'emploi : python use_check_etapes_et_termes.py nom_fig chemin_out directions  N_smooth [temps0,temps1]
#                                0                        1        2        3          4         5
# /!\ 3 : xyz ou xy ou x ou ... . all revient a xyz
# /!\ 4 : epaisseur du filtre pour adoucir les coubres 
#         (100 est ok pour le jdd de spectral_UO_enable_diphasique)
# /!\ 5 : [temps0,temps1] -> releve les donnees des outs entre les 2 temps
# /!\ 5 : [temps0,temps0] -> releve la donnee des outs a l'instant le plus proche de temps0
#         [-1,-1] -> releve tous les temps disponibles
########################################################################

message = "#############################################\n\
Tracons la QDM et les forces interfaciales au cours du temps \n\
##############################################"
print(message)

## Lecture des donnees utilisateur
nom_figure = sys.argv[1]
chemin_out = sys.argv[2]
directions = sys.argv[3]
N_smooth = float(sys.argv[4])
if directions == "all" : directions = "xyz"
t_deb, t_fin = float(sys.argv[5].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[5].rstrip("]").lstrip("[").split(",")[-1])
print("# t_deb = %s \n# t_fin = %s"%(str(t_deb),str(t_fin)))
########################################################################

## Creation des chemins, recuperation du temps
chemin_out = chemin_out+"*_check_etapes_et_termes.out"
print("shape")
print(numpy.loadtxt(glob.glob(chemin_out)[-1]).shape)
temps = numpy.loadtxt(glob.glob(chemin_out)[-1])[:,1]
########################################################################

# Travail sur les directions a considerer ##############################
liste_directions = []
if "x" in directions.lower() :
    liste_directions.append("x")
if "y" in directions.lower() :
    liste_directions.append("y")
if "z" in directions.lower() :
    liste_directions.append("z")


########################################################################

# Travail sur les iterations ###########################################
## Recuperation des iterations a considerer dans les fichiers out
# Si es instants renseignés ne sont pas inclus dans la plage de releves
#   -> l'instant initial est le prmeier instant de releve. L'instant final est le dernier instant de releve
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
########################################################################
# Pierre de Rosette pour savoir quelle colonne correspond à quelle grandeur
#              car dans le cas d'une reprise il se peut que l'on perde cette information
#              meme dans le cas d'un calcul initial, il a de fortes chances pour que ce
#              ne soit pas "bien" fait.
########################################################################
rosette = { "it":0,
            "t":1,
            "rk":2,
            "av_r_m_i" : {"x":1+2,  "y":1+3,  "z":1+4},   # rho_u_euler_av_rho_mu_ind
            "ap_r_m_i" : {"x":1+5,  "y":1+6,  "z":1+7},   # rho_u_euler_ap_rho_mu_ind
            "av_pred"  : {"x":1+8,  "y":1+9,  "z":1+10},  # rho_u_euler_av_prediction
            "ap_pred"  : {"x":1+11, "y":1+12, "z":1+13},  # rho_u_euler_ap_prediction
            "av_proj"  : {"x":1+14, "y":1+15, "z":1+16},  # rho_u_euler_av_projection
            "ap_proj"  : {"x":1+17, "y":1+18, "z":1+19},  # rho_u_euler_ap_projection
            "t_intf"   : {"x":1+20, "y":1+21, "z":1+22},  # terme_interfaces
            "t_conv"   : {"x":1+23, "y":1+24, "z":1+25},  # terme_convection
            "t_diff"   : {"x":1+26, "y":1+27, "z":1+28},  # terme_diffusion
            "t_pres"   : {"x":1+29, "y":1+30, "z":1+31}}  # terme_pression

color = {"x":"r","y":"g","z":"b"}
linestyle = ["solid", "dashed", "dotted", "dashdot"]
########################################################################

# Valeurs de references, propres au probleme
db = 1e-03
S = 4*numpy.pi*(db/2.)**2

## Lecture, selection, rangement et traitement des donnees
# chemin_out = "/volatile/RANDOM_TRIPERIO/a030/RE_ISOLE_143/RUN05/OUT/"

check = numpy.loadtxt(glob.glob(chemin_out)[-1])#[rosette["r_vel"]["x"]][it_deb:it_fin,1:]


## Sorties - not smoothed #############################################
fig_Pred, (ax_Pred) = plt.subplots(1,1,figsize=(10,10))
fig_Proj, (ax_Proj) = plt.subplots(1,1,figsize=(10,10))
fig_RMI, (ax_RMI) = plt.subplots(1,1,figsize=(10,10))
fig_PredProj, (ax_PredProj) = plt.subplots(1,1,figsize=(10,10))
fig_ProjRMI, (ax_ProjRMI) = plt.subplots(1,1,figsize=(10,10))
fig_all, (ax_all) = plt.subplots(1,1,figsize=(10,10))
fig_df, (ax_df) = plt.subplots(1,1,figsize=(10,10))
fig_pr, (ax_pr) = plt.subplots(1,1,figsize=(10,10))
fig_bilan, (ax_bilan) = plt.subplots(1,1,figsize=(10,10))


for i, direction in enumerate(liste_directions):
    ax_Pred.plot(temps,check[it_deb:it_fin+1,[rosette["av_pred"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" av. pred.")
    ax_Pred.plot(temps,check[it_deb:it_fin+1,[rosette["ap_pred"][direction]]], color=color[direction],linestyle=linestyle[1],label=direction+" ap. pred.")
                     
    ax_PredProj.plot(temps,check[it_deb:it_fin+1,[rosette["ap_pred"][direction]]], color=color[direction],linestyle=linestyle[0],label=direction+" ap. pred.")
    ax_PredProj.plot(temps,check[it_deb:it_fin+1,[rosette["av_proj"][direction]]], color=color[direction],linestyle=linestyle[1],label=direction+" av. proj.")
    # ax_PredProj.plot(temps,
                     # check[it_deb:it_fin+1,[rosette["ap_pred"][direction]]]-check[it_deb:it_fin+1,[rosette["av_proj"][direction]]]
                     # , color=color[direction],linestyle=linestyle[0],label=direction+" ap. pred. - av. proj.")

    ax_Proj.plot(temps,check[it_deb:it_fin+1,[rosette["av_proj"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" av. proj.")
    ax_Proj.plot(temps,check[it_deb:it_fin+1,[rosette["ap_proj"][direction]]], color=color[direction],linestyle=linestyle[1],label=direction+" ap. proj.")

    ax_ProjRMI.plot(temps,check[it_deb:it_fin+1,[rosette["ap_proj"][direction]]], color=color[direction],linestyle=linestyle[0],label=direction+" ap. proj.")
    ax_ProjRMI.plot(temps,check[it_deb:it_fin+1,[rosette["av_r_m_i"][direction]]], color=color[direction],linestyle=linestyle[1],label=direction+" av. rmi.")
    # ax_ProjRMI.plot(temps,
                     # check[it_deb:it_fin+1,[rosette["ap_proj"][direction]]]-check[it_deb:it_fin+1,[rosette["av_r_m_i"][direction]]]
                     # , color=color[direction],linestyle=linestyle[0],label=direction+" ap. pred. - av. proj.")

    ax_RMI.plot(temps,check[it_deb:it_fin+1,[rosette["av_r_m_i"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" av. rmi.")
    ax_RMI.plot(temps,check[it_deb:it_fin+1,[rosette["ap_r_m_i"][direction]]], color=color[direction],linestyle=linestyle[1],label=direction+" ap. rmi.")
    
    ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["av_pred"][direction]]], color="r"     ,linestyle="solid",label=direction+" av. pred.")
    ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["ap_pred"][direction]]], color="r",linestyle="dashed",label=direction+" ap. pred.")
    ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["av_proj"][direction]]], color="g"     ,linestyle="solid",label=direction+" av. proj.")
    ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["ap_proj"][direction]]], color="g",linestyle="dashed",label=direction+" ap. proj.")
    ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["av_r_m_i"][direction]]], color="b"     ,linestyle="solid",label=direction+" av. rmi.")
    ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["ap_r_m_i"][direction]]], color="b",linestyle="dashed",label=direction+" ap. rmi.")
    

ax_Pred.legend();ax_Pred.grid(True,'minor')
ax_Proj.legend();ax_Proj.grid(True,'minor')
ax_RMI.legend();ax_RMI.grid(True,'minor')
ax_PredProj.legend();ax_PredProj.grid(True,'minor')
ax_ProjRMI.legend();ax_ProjRMI.grid(True,'minor')
ax_all.legend();ax_all.grid(True,'minor')
ax_df.legend();ax_ProjRMI.grid(True,'minor')
ax_pr.legend();ax_ProjRMI.grid(True,'minor')
ax_bilan.legend();ax_bilan.grid(True,'minor')


ax_Pred.set_xlabel("temps (s)"    ,fontsize=22)
ax_Proj.set_xlabel("temps (s)"    ,fontsize=22)
ax_RMI.set_xlabel("temps (s)"     ,fontsize=22)
ax_PredProj.set_xlabel("temps (s)",fontsize=22)
ax_ProjRMI.set_xlabel("temps (s)" ,fontsize=22)
ax_all.set_xlabel("temps (s)"     ,fontsize=22)
ax_df.set_xlabel("temps (s)"      ,fontsize=22)
ax_pr.set_xlabel("temps (s)"      ,fontsize=22)
ax_bilan.set_xlabel("temps (s)"   ,fontsize=22)

ax_Pred.set_ylabel(r"$\rho^n u^{n}, \rho^n u^{*n}$ "      , fontsize=22)
ax_Proj.set_ylabel(r"$\rho^n u^{*n}, \rho^n u^{n+1}$ "    , fontsize=22)
ax_RMI.set_ylabel(r"$\rho^n u^{n+1}, \rho^{n+1} u^{n+1}$ ", fontsize=22)
ax_PredProj.set_ylabel(r" $\rho^n u^{*n}$ "               , fontsize=22)
ax_ProjRMI.set_ylabel(r"$\rho^n u^{n+1}$"                 , fontsize=22)
ax_all.set_ylabel(r"$\rho u$"                             , fontsize=22)
# ax_df.set_ylabel("Diffusion")
# ax_pr.set_ylabel("Pression")
# ax_bilan.set_ylabel("Bilan")

fig_Pred.savefig(nom_figure+"_Pred_"+directions+".png")
fig_Proj.savefig(nom_figure+"_Pred_"+directions+".png")
fig_RMI.savefig(nom_figure+"_RMI_"+directions+".png")
fig_PredProj.savefig(nom_figure+"_PredProj_"+directions+".png")
fig_ProjRMI.savefig(nom_figure+"_ProjRMI_"+directions+".png")
fig_all.savefig(nom_figure+"_all_"+directions+".png")
# fig_pr.savefig(nom_figure+"_pr_"+directions+".png")
# fig_df.savefig(nom_figure+"_df_"+directions+".png")
# fig_bilan.savefig(nom_figure+"_bilan_"+directions+".png")

