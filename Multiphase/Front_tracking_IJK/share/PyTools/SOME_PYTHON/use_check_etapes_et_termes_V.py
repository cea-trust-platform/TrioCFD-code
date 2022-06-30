# -*- coding: utf8

import os
import numpy
import subprocess
import scipy
from scipy import signal
# import pandas
import glob
import sys

import DNSTools3 as dtool
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
# matplotlib.use("Agg")
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

"""
    Construit les termes du bilan de qdm : convection; diffusion; interfaces; maj de rho; ...
    Travailler sur la fiche qui se trouve ici :
    /volatile/DERIVE_Y/BULLE_X/SIGMA_NON_NUL_RK3
    /volatile/EQUATION_ENERGIE/IJK_EXTENSION_3_DIPHASIQUE_CDP
"""

########################################################################
# mode d'emploi : python use_check_etapes_et_termes_V.py nom_fig chemin_out jdd directions  N_smooth [temps0,temps1]
#                                0                           1        2       3       4         5         6
# /!\ 3 : xyz ou xy ou x ou ... . all revient a xyz
# /!\ 4 : epaisseur du filtre pour adoucir les coubres 
#         (100 est ok pour le jdd de spectral_UO_enable_diphasique)
# /!\ 5 : [temps0,temps1] -> releve les donnees des outs entre les 2 temps
#         [temps0,temps0] -> releve la donnee des outs a l'instant le plus proche de temps0
#         [-1,-1] -> releve tous les temps disponibles
# exemple : python /volatile/RANDOM_TRIPERIO/SOME_PYTHON/use_check_etapes_et_termes_III.py fig_check OUT/ yz 15 [-1,-1]
########################################################################

message = "#############################################\n\
        Tracons les termes du bilan de QDM et les forces interfaciales au cours du temps \n\
        -> dt C\n\
        -> dt D\n\
        -> dt I\n\
        -> dt P\n\
        -> (r^{n+1} - r^{n}) u \n\
        ##############################################"
print(message)

## Lecture des donnees utilisateur
nom_figure = sys.argv[1]
chemin_out = sys.argv[2]
fichier_data = sys.argv[3]
directions = sys.argv[4]
N_smooth = float(sys.argv[5])
if directions == "all" : directions = "xyz"
t_deb, t_fin = float(sys.argv[6].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[6].rstrip("]").lstrip("[").split(",")[-1])
print("# t_deb = %s \n# t_fin = %s"%(str(t_deb),str(t_fin)))
########################################################################

## Creation des chemins, recuperation du temps
chemin_out = chemin_out+"*_check_etapes_et_termes.out"
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

########################################################################
# Pierre de Rosette pour savoir quelle colonne correspond à quelle grandeur
#              car dans le cas d'une reprise il se peut que l'on perde cette information
#              meme dans le cas d'un calcul initial, il a de fortes chances pour que ce
#              ne soit pas "bien" fait.
########################################################################
########################################################################
# Bonne note : 
# Dans le cadre de la quantification de l'effet de source et patch, on a s'intéresser aux grandeurs
#  o Convection : terme_moyen_convection_mass_solver_  =>  t_conv_mass_sol
#  o Diffusion  : terme_moyen_diffusion_mass_solver_   =>  t_diff_mass_sol
#  o Pression   : terme_pression_ter                   =>  t_pr_3
#  o Interfaces : terme_interfaces_af_mass_solver      =>  t_intf_af_ms
# qui sont homogenes a du/dt, tout comme les corrections
#  o Source     : qdm_source_                          =>  correction_qdm_colieaire_a_g 
#  o Patch      : qdm_patch_correction_                =>  correction_qdm_orthogonale_a_g
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

color = {"x":"r","y":"g","z":"b"}
linestyle = ["solid", "dashed", "dotted", "dashdot","loosely dashed"]
########################################################################
check = numpy.loadtxt(glob.glob(chemin_out)[-1])
########################################################################
if 2 in check[:,[rosette["rk"]]] :
    ## Le schema temporel de resolution est RK3
    print("# RK3 SIMULATION")
    rk3 = True
    indices_temps_rk3_2 = numpy.argwhere(check[:,[rosette["rk"]]]==2)[:,0]
    print("check[indices_temps_rk3_2[0],:3]",check[indices_temps_rk3_2[0],:3])
    temps_2 =  check[indices_temps_rk3_2,[rosette["t"]]]
    print("temps_2[0:5]",temps_2[0:5])
    # Version de MON_IJK anterieure au 24 juin : Le t en RK3 est mal ecrit
    # avec la mauvaise ecriture : 
    #     rk3 = 1 donne t1 = t + 5/12 dt
    #     rk3 = 0 donne t0 = t + 4/12 dt d'ou l'expression de dt
    mal_ecrit = False 
    if mal_ecrit:
        indices_temps_rk3_0 = numpy.argwhere(check[:,[rosette["rk"]]]==0)[:,0]
        indices_temps_rk3_1 = numpy.argwhere(check[:,[rosette["rk"]]]==1)[:,0]
        dt = 12. * (check[indices_temps_rk3_1 ,[rosette["t"]]] - check[indices_temps_rk3_0 ,[rosette["t"]]])
        # Puis les "bon" ti sont :
        #         t0 = t  + 4/12 dt        : inchange
        #         t1 = t0 + 5/12 dt        : change
        #         t2 = t1 + 3/12 dt = t+dt : change
        # /!\ La, on "peut" se contenter de faire comme ca car le out commence bien pas un rk_step = 0. Si ce n'etait pas le cas, 
        # /!\ il faudrai certainement se casser plus la tete pour bien faire. Dans la mesure ou c'est juste pour corriger une etourderie... on se casse pas plus 
        # Et paf, on remet les bonnes valeurs dans nos tableaux de temps, direct
        temps_0 = check[indices_temps_rk3_0 ,[rosette["t"]]]
        temps_1 = temps_0 + 5./12. * dt
        check[indices_temps_rk3_1 ,[rosette["t"]]] = temps_1
        temps_2 = temps_1 + 3./12. * dt
        check[indices_temps_rk3_2 ,[rosette["t"]]] = temps_2
        # temps   = temps_2

else:
    ## Le schema temporel de resolution est Euler explicite
    print("# EULER EXPLICIT SIMLATION")
    rk3 = False
temps = check[:,rosette["t"]]
########################################################################

########################################################################
# Travail sur les iterations ###########################################
# /!\ dans le cas de RK3, it_deb/fin sont plutot ligne_deb/fin
########################################################################
## Recuperation des iterations a considerer dans les fichiers out
# Si les instants renseignés ne sont pas inclus dans la plage de releves
#   -> l'instant initial est le premier instant de releve. L'instant final est le dernier instant de releve
if t_deb<=temps[0] : t_deb = temps[0]; print("Debut des releves : "+str(t_deb))
if (t_fin<=0 or t_fin>temps[-1]) : t_fin = temps[-1]; print("Fin des releves : "+str(t_fin))

# On associe les instants donnes a l'iteration correspondante
it_deb = numpy.argwhere(numpy.abs(temps-t_deb)==numpy.min(numpy.abs(temps-t_deb)))
#print("1 : ", it_deb)
if rk3:
    # On souhaite que les tableaux commencent par un rk_step = 2. Un point c'est out
    it=it_deb[0][0]
    vrai_it_deb = numpy.argwhere( numpy.abs(temps_2-temps[it]) == numpy.min(numpy.abs(temps_2-temps[it])))[0][0]
    #print("vrai_it_deb : ",vrai_it_deb)
    vrai_t_deb = temps_2[vrai_it_deb]
    #print("vrai_t_deb : ",vrai_t_deb)
    it_deb = numpy.argwhere(numpy.abs(temps-vrai_t_deb)==numpy.min(numpy.abs(temps-vrai_t_deb)))
#/!\ il faut sans doute faire la meme chose pour it_fin si le fichier out ne se termine pas par un rk = 2
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
print("# it_deb = %s \n# it_fin = %s"%(it_deb, it_fin))
temps = temps[it_deb:it_fin+1]
print("# t_deb = %s \n# t_fin = %s"%(temps[0], temps[-1]))
########################################################################


########################################################################
# On recupere uniquement les instants 'utiles' (en RK3, un pas sur trois avance réellement le temps)
# Dans IJK on a bien pris garde à faire un += quand il faut faire le += 
# /!\ ==> En schéma RK3 c'est la galère, on va construire les termes du bilan :
#                      terme_1 = (a0 * temre^{n}) + (a1 * terme^{n+1/3}) + (a2 * terme^{n+2/3})

if 2 in check[it_deb:it_fin,[rosette["rk"]]] :
    print("# RK3 SIMULATION")
    print("check[it_deb,:]",check[it_deb,:3])
    rk3 = True
    # On prend pour refernce de depart le rk = 2, il est attache a l'iteration precedent rk = 1
    indices_temps_rk3_2 = numpy.argwhere(check[it_deb:it_fin,[rosette["rk"]]]==2)[:,0] # je ne suis pas trop sur, mais ce doit etre une histoire que le decompte commence pas pareil pour tous   
    indices_temps_rk3_0 = numpy.argwhere(check[it_deb:it_fin,[rosette["rk"]]]==0)[:,0] # je ne suis pas trop sur, mais ce doit etre une histoire que le decompte commence pas pareil pour tous
    indices_temps_rk3_1 = numpy.argwhere(check[it_deb:it_fin,[rosette["rk"]]]==1)[:,0] # je ne suis pas trop sur, mais ce doit etre une histoire que le decompte commence pas pareil pour tous

    print("indices_temps_rk3_0[0:3]",indices_temps_rk3_0[0:3])
    print("indices_temps_rk3_1[0:3]",indices_temps_rk3_1[0:3])
    print("indices_temps_rk3_2[0:3]",indices_temps_rk3_2[0:3])

    print("indices_temps_rk3_0[-3:]",indices_temps_rk3_0[-3:])
    print("indices_temps_rk3_1[-3:]",indices_temps_rk3_1[-3:])
    print("indices_temps_rk3_2[-3:]",indices_temps_rk3_2[-3:])

    check_0 = check[indices_temps_rk3_0,...]
    check_1 = check[indices_temps_rk3_1,...]
    check_2 = check[indices_temps_rk3_2,...]
    
    print("check_0[0,0:3]",check_0[0,0:3])
    print("check_1[0,0:3]",check_1[0,0:3])
    print("check_2[0,0:3]",check_2[0,0:3])

    print("check_0.shape",check_0.shape)
    print("check_1.shape",check_1.shape)
    print("check_2.shape",check_2.shape)
    
    # Coefficients specifiques au RK3 de IJK, pour ecrire
    #     v = dt * (dv_0 * a_0 + dv_1 * a_1 + dv_2 * a_2)
    alpha_0 = 1./3. - 25./48. + (153*8.)/(9.*28*3)
    alpha_1 = 45/48. - (153*8.)/(128.*15)
    alpha_2 = 8/15.

    check_drho_u = (
            alpha_0 * (check_0[:,rosette["ru_av_rmi"]["x"]:rosette["ru_av_rmi"]["z"]+1] - check_0[:,rosette["ru_ap_rmi"]["x"]:rosette["ru_ap_rmi"]["z"]+1]) + 
            alpha_1 * (check_1[:,rosette["ru_av_rmi"]["x"]:rosette["ru_av_rmi"]["z"]+1] - check_1[:,rosette["ru_ap_rmi"]["x"]:rosette["ru_ap_rmi"]["z"]+1]) + 
            alpha_2 * (check_2[:,rosette["ru_av_rmi"]["x"]:rosette["ru_av_rmi"]["z"]+1] - check_2[:,rosette["ru_ap_rmi"]["x"]:rosette["ru_ap_rmi"]["z"]+1])  )
    print("check_drho_u.shape",check_drho_u.shape)

    # Constitution des sous pas de temps pour la resolution en Runge-Kutta (1/3, 5/12, 1/4)
    # Le premier rk_step vaut 2. Le dernier rk_step vaut 2.
    # Finalement, on peut s'en passer des sous_dt. Gardons-les desfois qu'il me les faudrai
    sous_dt_0 = (check_0[0:,[rosette["t"]]] - check_2[0:,[rosette["t"]]]).flatten()
    sous_dt_1 = (check_1[0:,[rosette["t"]]] - check_0[0:,[rosette["t"]]]).flatten()
    sous_dt_2 = (check_2[0:,[rosette["t"]]] - check_1[0:,[rosette["t"]]]).flatten()


    print("sous_dt_0[0:5]",sous_dt_0[0:5])
    print("sous_dt_1[0:5]",sous_dt_1[0:5])
    print("sous_dt_2[0:5]",sous_dt_2[0:5])
    # temps = check_0[:,[rosette["t"]]]
else :
    print("# EULER EXPLICIT SIMLATION")
    rk3 = False
    check_0 = check[it_deb:it_fin,...]
    temps = check_0[it_deb:it_fin,[rosette["t"]]]
print("# Again the time : ")
print("# t_deb = %s \n# t_fin = %s"%(temps[0], temps[-1]))

# Valeurs de references, propres au probleme
db = 1e-03
S = 4*numpy.pi*(db/2.)**2
rl = dtool.getParam(fichier_data, "rho_liquide")
rv = dtool.getParam(fichier_data, "rho_vapeur")

## Lecture, selection, rangement et traitement des donnees
# chemin_out = "/volatile/RANDOM_TRIPERIO/a030/RE_ISOLE_143/RUN05/OUT/"


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
# Declaration figures
fig_t_intf_bis, (ax_t_intf_bis) = plt.subplots(1,1,figsize=(10,10))
fig_t_intf, (ax_t_intf) = plt.subplots(1,1,figsize=(10,10))
fig_t_conv, (ax_t_conv) = plt.subplots(1,1,figsize=(10,10))
fig_t_diff, (ax_t_diff) = plt.subplots(1,1,figsize=(10,10))
fig_t_rho,  (ax_t_rho)  = plt.subplots(1,1,figsize=(10,10)) # Effet du terme (rho^{n+1} - rho^{n}) * u^{n+1}
fig_t_pr_2, (ax_t_pr_2) = plt.subplots(1,1,figsize=(10,10)) # grad(p)
fig_t_pr_3, (ax_t_pr_3) = plt.subplots(1,1,figsize=(10,10)) # 1/\rho grad(p)
fig_t_pr_d, (ax_t_pr_d) = plt.subplots(1,1,figsize=(10,10))
fig_all,    (ax_all)    = plt.subplots(1,1,figsize=(10,10))
fig_ru,     (ax_ru)     = plt.subplots(1,1,figsize=(10,10))
fig_bilan,  (ax_bilan)  = plt.subplots(1,1,figsize=(10,10))

fig_t_intf_ms, (ax_t_intf_ms) = plt.subplots(1,1,figsize=(10,10))
fig_t_conv_ms, (ax_t_conv_ms) = plt.subplots(1,1,figsize=(10,10))
fig_t_diff_ms, (ax_t_diff_ms) = plt.subplots(1,1,figsize=(10,10))
fig_all_ms,    (ax_all_ms)    = plt.subplots(1,1,figsize=(10,10))

fig_t_intf.suptitle(r"$ dt \times \langle \frac{\mathcal{F}_i}{\mathcal{V} }\rangle$", fontsize=24)
fig_t_intf_ms.suptitle(r"$ dt \times \langle \frac{\mathcal{F}_i}{\tilde{\rho} \mathcal{V} } \rangle$", fontsize=24)

# Titres figures
option_convection = "non_conservative_simple"
if option_convection == "non_conservative_simple":
    fig_t_conv.suptitle(r"$ dt \times \langle \tilde{\rho} \mathcal{S} \left[ -u_i u_j \right] \rangle \frac{1}{\mathcal{V}}$", fontsize=24)
    fig_t_conv_ms.suptitle(r"$ dt \times \frac{ \langle \tilde{\rho} \mathcal{S} \left[ -u_i u_j \right] \rangle }{\tilde{\rho}} \frac{1}{\mathcal{V}}$", fontsize=24)
elif option_convection == "non_conservative_rhou":
    fig_t_conv.suptitle(r"$ dt \times \langle \mathcal{S} \left[ -u_i \rho uj \right] + u_i \partial_j \rho u_j  \rangle \frac{1}{\mathcal{V}}$", fontsize=24)
elif option_convection == "conservative":
    fig_t_conv.suptitle(r"$ dt \times \langle \mathcal{S} \left[ -u_i \rho uj \right] \rangle \frac{1}{\mathcal{V}}$", fontsize=24)

fig_t_diff.suptitle(r"$ dt \times \langle \mathcal{S} \left[ \mu s_{ij} \right] + \mathcal{C} \rangle \frac{1}{\mathcal{V}} - \langle \mathcal{C} \rangle \frac{1}{\mathcal{V}}$", fontsize=24)
fig_t_diff_ms.suptitle(r"$ dt \times \frac{ \langle \mathcal{S} \left[ \mu s_{ij} \right] + \mathcal{C} \rangle \frac{1}{\mathcal{V}} - \langle \mathcal{C} \rangle \frac{1}{\mathcal{V}} }{\tilde{\rho}}$", fontsize=24)
u=0; rho_u = 1
if u:
    fig_t_pr_2.suptitle(r"$ dt \times \nabla{p}$", fontsize=24)
    fig_t_pr_3.suptitle(r"$ dt \times \frac{1}{\rho} \nabla{p}$", fontsize=24)
    fig_t_pr_d.suptitle(r"$ \tilde{\rho}^* dt \times \nabla{p} _{1/2} - dt \times \frac{1}{\rho} \nabla{p}$", fontsize=24)
elif rho_u:
    fig_t_pr_2.suptitle(r"$ dt \times \nabla{p}$", fontsize=24)
    fig_t_pr_3.suptitle(r"$ dt \times \nabla{p}$, true", fontsize=24)
    fig_t_pr_d.suptitle(r"$ dt \times \nabla{p}$, diff", fontsize=24)
fig_t_rho.suptitle(r"$\rho^{n} u^{n+1} -  \rho^{n+1} u^{n+1}$", fontsize=24)

fig_all.suptitle(r"Somme des termes", fontsize=24)
fig_all_ms.suptitle(r"Somme des termes, after m.s.", fontsize=24)
fig_ru.suptitle(r"Detail des termes", fontsize=24)


# Construction du pas de temps de simulation
dt = check_0[1:,[rosette["t"]]] - check_0[:-1,[rosette["t"]]]
dt = numpy.insert(dt,0,check_0[1,[rosette["t"]]])

# On reforme le temps : les sous-pas de temps n'existent plus. Ils ont existé uniquement pour re-construire les termes totaux 
if rk3:
    it_deb_sub = vrai_it_deb # defini la ou on se debat avec les it_deb et tout, au debut
    it_fin_sub = numpy.argwhere( numpy.abs(temps_2-t_fin == numpy.min(numpy.abs(temps_2-t_fin))))[0][0]
    temps_plot = temps_2[it_deb_sub:it_fin_sub]
else:
    temps_plot = temps

for i, direction in enumerate(liste_directions):
    # Construction des termes du bilan. Homogenes a du/dt pour les quatre premiers
    # Le dernier est homogene a du. Il a un traitement (et une definition) qui dépend du schema temporel de resolution
    intf = check_0[:,[rosette["t_intf_bf_ms"][direction]]].flatten()
    conv = check_0[:,[rosette["t_conv"][direction]]].flatten()
    pr_2 = check_0[:,[rosette["t_pr_2"][direction]]].flatten()
    pr_3 = check_0[:,[rosette["t_pr_3"][direction]]].flatten()
    # On s'attedns a avoir pr_d = 0 partout, tout le temps
    pr_d = (1./rl + 1./rv)/2.*pr_2 - pr_3
    diff = check_0[:,[rosette["t_diff"][direction]]].flatten()
    drho_u = (check_0[:,[rosette["ru_av_rmi"][direction]]] - check_0[:,[rosette["ru_ap_rmi"][direction]]]).flatten()
    d_rhou = (check_0[1:,[rosette["ru_ap_rmi"][direction]]] - check_0[:-1,[rosette["ru_ap_rmi"][direction]]])  
    d_rhou = numpy.insert(d_rhou,0,0)
    # Grandeurs apres etre passees au mass solver
    intf_ms = check_0[:,[rosette["t_intf_af_ms"][direction]]].flatten()
    conv_ms = check_0[:,[rosette["t_conv_mass_sol"][direction]]].flatten()
    diff_ms = check_0[:,[rosette["t_diff_mass_sol"][direction]]].flatten()



    if rk3 :
        intf_1, intf_2 = check_1[:,[rosette["t_intf_bf_ms"][direction]]].flatten(), check_2[:,[rosette["t_intf_bf_ms"][direction]]].flatten()
        conv_1, conv_2 = check_1[:,[rosette["t_conv"][direction]]].flatten()      , check_2[:,[rosette["t_conv"][direction]]].flatten()      
        pr_2_1, pr_2_2 = check_1[:,[rosette["t_pr_2"][direction]]].flatten()      , check_2[:,[rosette["t_pr_2"][direction]]].flatten()      
        pr_3_1, pr_3_2 = check_1[:,[rosette["t_pr_3"][direction]]].flatten()      , check_2[:,[rosette["t_pr_3"][direction]]].flatten()      
        pr_d_1, pr_d_2 = pr_2_1 - pr_2_1                                          , pr_2_2 - pr_3_2     
        diff_1, diff_2 = check_1[:,[rosette["t_diff"][direction]]].flatten()      , check_2[:,[rosette["t_diff"][direction]]].flatten()      
        # Grandeurs apres etre passees au mass solver
        intf_ms_1, intf_ms_2 = check_1[:,[rosette["t_intf_af_ms"][direction]]].flatten()      , check_2[:,[rosette["t_intf_af_ms"][direction]]].flatten()      
        conv_ms_1, conv_ms_2 = check_1[:,[rosette["t_conv_mass_sol"][direction]]].flatten()      , check_2[:,[rosette["t_conv_mass_sol"][direction]]].flatten()      
        diff_ms_1, diff_ms_2 = check_1[:,[rosette["t_diff_mass_sol"][direction]]].flatten()      , check_2[:,[rosette["t_diff_mass_sol"][direction]]].flatten()      

        drho_u = check_drho_u[:,i] # (check[:,[rosette["ru_av_rmi"][direction]]] - check[:,[rosette["ru_ap_rmi"][direction]]]).flatten()
        drho_u_0 =  (check_0[1:,[rosette["ru_ap_rmi"][direction]]] - check_0[:-1,[rosette["ru_ap_rmi"][direction]]]) 
        # drho_u = numpy.insert(drho_u,0,0)    
        print("drho_u.shape",drho_u.shape)
        print("pressions",pr_3,pr_3_1,pr_3_2)
        print("conv*sous_dt_0.shape",(conv*sous_dt_0).shape)
        print("conv.shape,sous_dt_0.shape",conv.shape,sous_dt_0.shape)
        intf = intf*alpha_0 + intf_1*alpha_1 + intf_2*alpha_2
        conv = conv*alpha_0 + conv_1*alpha_1 + conv_2*alpha_2
        pr_2 = pr_2*alpha_0 + pr_2_1*alpha_1 + pr_2_2*alpha_2
        pr_3 = pr_3*alpha_0 + pr_3_1*alpha_1 + pr_3_2*alpha_2
        pr_d = pr_d*alpha_0 + pr_d_1*alpha_1 + pr_d_2*alpha_2
        diff = diff*alpha_0 + diff_1*alpha_1 + diff_2*alpha_2
        # Grandeurs apres etre passe au mass solver
        intf_ms = intf_ms*alpha_0 + intf_ms_1*alpha_1 + intf_ms_2*alpha_2
        conv_ms = conv_ms*alpha_0 + conv_ms_1*alpha_1 + conv_ms_2*alpha_2
        diff_ms = diff_ms*alpha_0 + diff_ms_1*alpha_1 + diff_ms_2*alpha_2
    
        pr_d_dt = pr_d*dt
    # else:
    intf_dt = intf*dt
    conv_dt = conv*dt
    pr_2_dt = pr_2*dt
    pr_3_dt = pr_3*dt
    diff_dt = diff*dt
    # Grandeurs apres etre passe au mass solver
    intf_ms_dt = intf_ms*dt
    conv_ms_dt = conv_ms*dt
    diff_ms_dt = diff_ms*dt
    # print("intf_ms_dt",intf_ms_dt)
    ###################################################
    # Construction du bilan, homogene a rho du pour les variables "normales" et homognee a du pour les variables en _ms (: mass solver)
    # Etre exactement fidele au code dans la reconstruction demande trop de manipulations
    # Pour le bilan on prend donc la grandeur grad(p)
    # Si on veut etre tres fidele a IJK, il faut calculer t_conv, t_diff, t_inf APRES QU'IL SOIENT PASSÉS DANS LE MASS SOLVER
    # pour ce faire, il faudrai passer {t_conv} au mass solver, puis {t_diff + t_conv} au mass solveur. On en deduit {t_diff} apres mass solver
    # et ainsi de suite pour les autres termes qui sont ajoutes AVANT la contribution de la gravite
    # ATTENTION les interfaces sont ajoutees directement a la velocity.
    # Dans ce cas, il va falloir creer de nouveaux champs...
    
    # bilan = (intf+conv+pr_2+diff)*dt + drho_u
    bilan = intf_dt+conv_dt+pr_2_dt+diff_dt + drho_u                          # @ d(\rho u)/dt
    bilan_ms = intf_ms_dt+conv_ms_dt+pr_3_dt+diff_ms_dt + drho_u*(2./(rl+rv))  # @ du/dt   (apres mass solver)
    # terme interfaces
    ax_t_intf_bis.plot(temps_plot,intf_dt, color=color[direction]     ,linestyle=linestyle[3],label=direction+" interfaces bf. ms : I")
    ax_t_intf_ms.plot(temps_plot,intf_ms_dt, color=color[direction]     ,linestyle=linestyle[1],label=direction+" interfaces af. ms")

    # terme convection
    ax_t_conv.plot(temps_plot,conv_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" convection")
    ax_t_conv_ms.plot(temps_plot,conv_ms_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" conv. after m.s.")

    # terme pression
    ax_t_pr_2.plot(temps_plot,pr_2_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" pression")
    ax_t_pr_3.plot(temps_plot,pr_3_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" pression")
    ax_t_pr_d.plot(temps_plot,pr_d_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" pression")

    # terme diffusion
    ax_t_diff.plot(temps_plot,diff_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" diffusion")
    ax_t_diff_ms.plot(temps_plot,diff_ms_dt, color=color[direction]     ,linestyle=linestyle[0],label=direction+" diff. after m.s.")

    # terme de transport des interfaces
    ax_t_rho.plot(temps_plot,drho_u, color=color[direction]     ,linestyle=linestyle[0],label=direction+r" $ (d_n \rho) \times u^{n+1}$")

    # Bilan : convection et diffusion sur dv lors de la prediction, interfaces sur velocity lors de la prediciton, pression sur velocity lors de la projection.
    # TODO : Ajouter terme_gravite dans IJK de la mm facon que les autres termes : dvelocity - convection - diffusion    
    ax_all.plot(temps_plot,bilan, color=color[direction]     ,linestyle="solid",label=direction+r" $(\mathcal{I} + \mathcal{C} + \mathcal{P} + \mathcal{D})\times dt + (\partial \rho) u$")
    ax_all_ms.plot(temps_plot,bilan_ms, color=color[direction]     ,linestyle="solid",label=direction+r" $(\mathcal{I} + \mathcal{C} + \mathcal{P} + \mathcal{D})\times dt + (\partial \rho) u$, after m.s.")
    print("temps_plot.shape, d_rhou.shape",temps_plot.shape, d_rhou.shape)
    ax_all.plot(temps_plot, d_rhou, color=color[direction]     ,linestyle="dashed",label=direction+r"$\partial (\rho u) $")

    # Detail des termes
    ax_ru.plot(temps_plot, conv_dt ,color=color[direction],linestyle=linestyle[0],label=direction+" convection")
    ax_ru.plot(temps_plot, intf_dt ,color=color[direction],linestyle=linestyle[1],label=direction+" interfaces")
    ax_ru.plot(temps_plot, pr_2_dt ,color=color[direction],linestyle=linestyle[2],label=direction+" pression") 
    ax_ru.plot(temps_plot, pr_3_dt ,color=color[direction],linestyle=linestyle[2],alpha=0.5,linewidth=3,label=direction+" pression") 
    ax_ru.plot(temps_plot, diff_dt ,color=color[direction],linestyle=linestyle[3],label=direction+" diffusion") 
    ax_ru.scatter(temps_plot, drho_u ,color=color[direction],s=5,label=direction+" transport chi") 
    # ax_ru.scatter(temps_plot, d_rhou ,color=color[direction],s = 7,label=direction+r"$ \partial_t (\rho u)$") 

ax_t_intf_bis.legend();ax_t_intf_bis.grid(True,'minor')
ax_t_intf.legend();ax_t_intf.grid(True,'minor')
ax_t_intf_ms.legend();ax_t_intf_ms.grid(True,'minor')
ax_t_conv.legend();ax_t_conv.grid(True,'minor')
ax_t_conv_ms.legend();ax_t_conv_ms.grid(True,'minor')
ax_t_diff.legend();ax_t_diff.grid(True,'minor')
ax_t_diff_ms.legend();ax_t_diff_ms.grid(True,'minor')
ax_all.legend();ax_all.grid(True,'minor')
ax_all_ms.legend();ax_all_ms.grid(True,'minor')
ax_ru.legend();ax_ru.grid(True,'minor')


ax_t_intf_bis.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_intf.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_intf_ms.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_conv.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_conv_ms.set_xlabel("temps (s)"    ,fontsize=22)
ax_t_diff.set_xlabel("temps (s)"     ,fontsize=22)
ax_t_diff_ms.set_xlabel("temps (s)"     ,fontsize=22)
ax_all.set_xlabel("temps (s)"     ,fontsize=22)
ax_all_ms.set_xlabel("temps (s)"     ,fontsize=22)
ax_t_pr_2.set_xlabel("temps (s)"      ,fontsize=22)
ax_t_pr_3.set_xlabel("temps (s)"      ,fontsize=22)
ax_t_pr_d.set_xlabel("temps (s)"      ,fontsize=22)
ax_ru.set_xlabel("temps (s)"   ,fontsize=22)

ax_t_intf_bis.set_ylabel(r"$\int \sigma^{n} \kappa^{n} n $ "      , fontsize=22)
ax_t_intf.set_ylabel(r"$\int \sigma^{n} \kappa^{n} n $ "      , fontsize=22)
ax_t_intf_ms.set_ylabel(r"$\frac{1}{\tilde{\rho}} \int \sigma^{n} \kappa^{n} n $ "      , fontsize=22)
ax_t_conv.set_ylabel(r"$u^{*n} \nabla u^{*n}$ "    , fontsize=22)
ax_t_conv_ms.set_ylabel(r"$ \frac{1}{\tilde{\rho}} u^{*n} \nabla u^{*n}$ "    , fontsize=22)
ax_t_diff.set_ylabel(r"$2 \nu \nabla s_{ij}$ ", fontsize=22)
ax_t_diff_ms.set_ylabel(r"$\frac{1}{\tilde{\rho}} 2 \nu \nabla s_{ij}$ ", fontsize=22)
ax_t_pr_2.set_ylabel(r"$\frac{1}{\rho} \nabla p$ ", fontsize=22)
ax_t_pr_3.set_ylabel(r"$\frac{1}{\rho} \nabla p$ true ", fontsize=22)
ax_t_pr_d.set_ylabel(r"$\frac{1}{\rho} \nabla p$ diff ", fontsize=22)
ax_all.set_ylabel(r"$\rho u$"                             , fontsize=22)
ax_all.set_ylabel(r"$ u$, bilan"                             , fontsize=22)
ax_ru.set_ylabel(r"detail des termes"                             , fontsize=22)

# plt.show()
fig_t_intf_bis.savefig(nom_figure+"_t_intf_bis_"+directions+".png")
# fig_t_intf.savefig(nom_figure+"_t_intf_"+directions+".png")
fig_t_intf_ms.savefig(nom_figure+"_t_intf_ms"+directions+".png")
fig_t_conv.savefig(nom_figure+"_t_conv_"+directions+".png")
fig_t_conv_ms.savefig(nom_figure+"_t_conv_ms"+directions+".png")
fig_t_diff.savefig(nom_figure+"_t_diff_"+directions+".png")
fig_t_diff_ms.savefig(nom_figure+"_t_diff_ms"+directions+".png")
fig_t_rho.savefig(nom_figure+"_t_rho_"+directions+".png")
fig_t_pr_2.savefig(nom_figure+"_t_pr_2_"+directions+".png")
fig_t_pr_3.savefig(nom_figure+"_t_pr_3_"+directions+".png")
fig_t_pr_d.savefig(nom_figure+"_t_pr_d_"+directions+".png")
fig_all.savefig(nom_figure+"_all_"+directions+".png")
fig_all_ms.savefig(nom_figure+"_all_ms_"+directions+".png")
fig_ru.savefig(nom_figure+"_detail_"+directions+".png")

