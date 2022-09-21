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
            "t_pres"   : {"x":1+29, "y":1+30, "z":1+31},  # terme_pression
            "t_intf_bf_ms_bis"   : {"x":33, "y":34, "z":35},  # terme_interfaces_bf_mass_solveur bis
            "t_intf_bf_ms"       : {"x":36, "y":37, "z":38},  # terme_interfaces_bf_mass_solveur
            "t_intf_af_ms"       : {"x":39, "y":40, "z":41}}  # terme_interfaces_af_mass_solveur

color = {"x":"r","y":"g","z":"b"}
linestyle = ["solid", "dashed", "dotted", "dashdot"]
########################################################################

# Valeurs de references, propres au probleme
db = 1e-03
S = 4*numpy.pi*(db/2.)**2

## Lecture, selection, rangement et traitement des donnees
# chemin_out = "/volatile/RANDOM_TRIPERIO/a030/RE_ISOLE_143/RUN05/OUT/"

check = numpy.loadtxt(glob.glob(chemin_out)[-1])#[rosette["r_vel"]["x"]][it_deb:it_fin,1:]
########################################################################
# ////// Extraits d'IJK //////
# Termes calcules pendant l'étape de prédiction 

# convection (a partir de dv)
# terme_convection[i] = calculer_v_moyen(d_velocity_[i])/volume_cell_uniforme;

# diffusion (a partir de dv et du terme de convection)
# terme_diffusion[i] = calculer_v_moyen(d_velocity_[i])/volume_cell_uniforme - terme_convection[i];

# interfaces (directement le terme calculé)
# terme_interfaces_af_mass_solver[dir] = calculer_v_moyen(terme_source_interfaces_ns_[dir]);
# terme_interfaces_bf_mass_solver[dir] = calculer_v_moyen(terme_source_interfaces_ns_[dir]);

# interfaces (a partir de dv et des autres termes)
# terme_interfaces_bf_mass_solver_bis[dir] = calculer_v_moyen(d_velocity_[dir])/volume_cell_uniforme;

## Sorties - not smoothed #############################################
fig_t_intf_bis, (ax_t_intf_bis) = plt.subplots(1,1,figsize=(10,10))
fig_t_intf, (ax_t_intf) = plt.subplots(1,1,figsize=(10,10))
fig_t_conv, (ax_t_conv) = plt.subplots(1,1,figsize=(10,10))
fig_t_diff, (ax_t_diff) = plt.subplots(1,1,figsize=(10,10))
fig_t_pres, (ax_t_pres) = plt.subplots(1,1,figsize=(10,10))
fig_all, (ax_all) = plt.subplots(1,1,figsize=(10,10))
fig_t_convt_diff, (ax_t_convt_diff) = plt.subplots(1,1,figsize=(10,10))
fig_t_intft_conv, (ax_t_intft_conv) = plt.subplots(1,1,figsize=(10,10))
fig_pr, (ax_pr) = plt.subplots(1,1,figsize=(10,10))
fig_bilan, (ax_bilan) = plt.subplots(1,1,figsize=(10,10))


for i, direction in enumerate(liste_directions):
    ax_t_intf_bis.plot(temps,check[it_deb:it_fin+1,[rosette["t_intf_bf_ms_bis"][direction]]], color="orange"     ,linestyle=linestyle[0],label=direction+" interfaces bf. ms")
    
    ax_t_intf.plot(temps,check[it_deb:it_fin+1,[rosette["t_intf_bf_ms"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" interfaces bf. ms")
    ax_t_intf.plot(temps,check[it_deb:it_fin+1,[rosette["t_intf_af_ms"][direction]]], color=color[direction]     ,linestyle=linestyle[1],label=direction+" interfaces af. ms")
    ax_t_intf.scatter(temps,check[it_deb:it_fin+1,[rosette["t_intf_af_ms"][direction]]], color="cyan"     ,label=direction+" interfaces af. ms")
    
    ax_t_conv.plot(temps,check[it_deb:it_fin+1,[rosette["t_conv"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" convection")
    
    ax_t_pres.plot(temps,check[it_deb:it_fin+1,[rosette["t_pres"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" pression")
    
    ax_t_diff.plot(temps,check[it_deb:it_fin+1,[rosette["t_diff"][direction]]], color=color[direction]     ,linestyle=linestyle[0],label=direction+" diffusion")
    
    # Bilan, à l'issu de la prediction, mais sans le "vrai" terme de pression calculé (il est determine lors de la projection)
    intf = check[it_deb:it_fin+1,[rosette["t_intf_af_ms"][direction]]]
    conv = check[it_deb:it_fin+1,[rosette["t_conv"][direction]]]
    pres = check[it_deb:it_fin+1,[rosette["t_pres"][direction]]]
    diff = check[it_deb:it_fin+1,[rosette["t_diff"][direction]]]
    ax_all.plot(temps,intf+conv+pres+diff, color="r"     ,linestyle="solid",label=direction+" sum")
    # ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["ap_pred"][direction]]], color="r",linestyle="dashed",label=direction+" ap. pred.")
    # ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["av_proj"][direction]]], color="g"     ,linestyle="solid",label=direction+" av. proj.")
    # ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["ap_proj"][direction]]], color="g",linestyle="dashed",label=direction+" ap. proj.")
    # ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["av_r_m_i"][direction]]], color="b"     ,linestyle="solid",label=direction+" av. rmi.")
    # ax_all.plot(temps,check[it_deb:it_fin+1,[rosette["ap_r_m_i"][direction]]], color="b",linestyle="dashed",label=direction+" ap. rmi.")
    

ax_t_intf_bis.legend();ax_t_intf_bis.grid(True,'minor')
ax_t_intf.legend();ax_t_intf.grid(True,'minor')
ax_t_conv.legend();ax_t_conv.grid(True,'minor')
ax_t_diff.legend();ax_t_diff.grid(True,'minor')
ax_t_intft_conv.legend();ax_t_intft_conv.grid(True,'minor')
ax_t_convt_diff.legend();ax_t_convt_diff.grid(True,'minor')
ax_all.legend();ax_all.grid(True,'minor')
ax_t_pres.legend();ax_t_convt_diff.grid(True,'minor')
ax_pr.legend();ax_t_convt_diff.grid(True,'minor')
ax_bilan.legend();ax_bilan.grid(True,'minor')


ax_t_intf_bis.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_intf.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_conv.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_diff.set_xlabel("temps (s)"     ,fontsize=22)
ax_t_intft_conv.set_xlabel("temps (s)",fontsize=22)
ax_t_convt_diff.set_xlabel("temps (s)" ,fontsize=22)
ax_all.set_xlabel("temps (s)"     ,fontsize=22)
ax_t_pres.set_xlabel("temps (s)"      ,fontsize=22)
ax_pr.set_xlabel("temps (s)"      ,fontsize=22)
ax_bilan.set_xlabel("temps (s)"   ,fontsize=22)

ax_t_intf_bis.set_ylabel(r"$\int \sigma^{n} \kappa^{n} n $ "      , fontsize=22)
ax_t_intf.set_ylabel(r"$\int \sigma^{n} \kappa^{n} n $ "      , fontsize=22)
ax_t_conv.set_ylabel(r"$u^{*n} \nabla u^{*n}$ "    , fontsize=22)
ax_t_diff.set_ylabel(r"$2 \nu \nabla s_{ij}$ ", fontsize=22)
ax_t_pres.set_ylabel(r"$\frac{1}{\rho} \nabla p$ ", fontsize=22)
ax_all.set_ylabel(r"$\rho u$"                             , fontsize=22)
# ax_t_intft_conv.set_ylabel(r" $\rho^n u^{*n}$ "               , fontsize=22)
# ax_t_convt_diff.set_ylabel(r"$\rho^n u^{n+1}$"                 , fontsize=22)
# ax_t_pres.set_ylabel("Diffusion")
# ax_pr.set_ylabel("Pression")
# ax_bilan.set_ylabel("Bilan")

fig_t_intf_bis.savefig(nom_figure+"_t_intf_bis_"+directions+".png")
fig_t_intf.savefig(nom_figure+"_t_intf_"+directions+".png")
fig_t_conv.savefig(nom_figure+"_t_conv_"+directions+".png")
fig_t_diff.savefig(nom_figure+"_t_diff_"+directions+".png")
fig_t_pres.savefig(nom_figure+"_t_pres_"+directions+".png")
fig_all.savefig(nom_figure+"_all_"+directions+".png")
# fig_pr.savefig(nom_figure+"_pr_"+directions+".png")
# fig_bilan.savefig(nom_figure+"_bilan_"+directions+".png")
# fig_t_intft_conv.savefig(nom_figure+"_Predt_conv_"+directions+".png")
# fig_t_convt_diff.savefig(nom_figure+"_t_convt_diff_"+directions+".png")

