# -*- coding: utf8

import os
import numpy as np
import subprocess
import scipy
from scipy import signal
import pandas
import glob
import sys

import scipy.stats

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist

from depSondes_classes_et_fonctions import *
import DNSTools3 as dtool

from scipy.integrate import simps

import math as mt

"""
    Correlations de la grandeur-sonde demandé
"""

##############################################################################
# mode d'emploi : python correlation_vIII.py chemin_son  pattern  quelle_grandeur  N_smooth
#                            0                1            2       3                 4   
#                        [direction_sondes]   [direction_force]  [temps0,temps1]
#                               5                  6                   7
# /!\ 1 : chemin_son -> Il FAUT que le chemin finisse par sont "/"
# /!\ 2 : pattern -> Pour traiter moins de sondes, c'est soit **, *XX* ou XX sont les
#         correspondant a une sonde
# /!\ 3 : Marqueur de la grandeur pour laquelle trcer les correlations.
# /!\ 4 : Epaisseur du filtre pour moyenner le signal temporel 
#         par exemple, N_smooth = 4 est le plus petit N_smooth qui donne un résultat visuellement ok
#         pour le fiche /volatile/FFTW/THI_FFT/Rapport_validation/spectral_UO
# /!\ 5 : direction_sondes -> DOIT etre une liste, meme si on ne traite que les sondes X
# /!\ 5 : direction_forces -> DOIT etre une liste, meme si on ne traite que les vitesses X
# /!\ 7 : [temps0,temps1] -> releve les donnees des outs entre les 2 temps
# /!\ 7 : [temps0,temps0] -> releve la donnee des outs a l'instant le plus proche de temps0
#         [-1,-1] -> releve tous les temps disponibles
# Exemple : 
# $ python correlation_vIII.py SON/  '**'  V  15 [X,Y,Z] [X,Y,Z] [-1,-1]
##############################################################################

message = "###########################################################\n\
# Evolution des autocorrelations temporelles et spatiales du forcage impose\n\
###"

"""
Etape 0 : lecture paramètres utilisateur, traitement de ces parametres
           initialisation des parametres fxes du probleme (lisibles dans le jdd en vrai)
           instruction des choses a faire 'do_the_[bidule]'
"""
## Lecture des donnees utilisateur
print("#######")
print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(sys.argv[5])
print(sys.argv[6])
print(sys.argv[7])
print("#######")

chemin_son = sys.argv[1]
pattern = sys.argv[2]
quelle_grandeur = sys.argv[3]
N_smooth = int(sys.argv[4])
direction_sondes = sys.argv[5].rstrip("]").lstrip("[").split(",") 
# Ajustement des sondes a traiter 
for i, ds in enumerate(direction_sondes):
    direction_sondes[i] += pattern

print("# Sondes traittees : \n# {0}".format(direction_sondes))
direction_forces = sys.argv[6].rstrip("]").lstrip("[").split(",")
# Ajustement des grandeurs a traiter
for i, dv in enumerate(direction_forces):
    direction_forces[i] = quelle_grandeur+dv 

print("# Vitesses traitees : \n# {0}".format(direction_forces))
t_deb, t_fin = float(sys.argv[7].rstrip("]").lstrip("[").split(",")[-2]),float(sys.argv[7].rstrip("]").lstrip("[").split(",")[-1])
print("# Demande :\n# t_deb = %s \n# t_fin = %s"%(str(t_deb),str(t_fin)))
########################################################################

## Creation des chemins, recuperation du temps
sonde00 = chemin_son+'*_S_'+direction_sondes[0]+'_'+direction_forces[0]+'.son'
print("sonde00",sonde00)
temps = np.loadtxt(glob.glob(sonde00)[-1])[5:,0]
########################################################################

## Recuperation des iterations a considerer dans les fichiers out
# Si es instants renseignés ne sont pas inclus dans la plage de releves
#   -> l'instant initial est le prmeier instant de releve. L'instant final est le dernier instant de releve
print("# Realise :")
if (t_deb<=temps[0]) : t_deb = temps[0]
if (t_fin<=0 or t_fin>temps[-1]) : t_fin = temps[-1]
# On associe les instants donnes a l'iteration correspondante
it_deb = np.argwhere(np.abs(temps-t_deb)==np.min(np.abs(temps-t_deb)))
it_fin = np.argwhere(np.abs(temps-t_fin)==np.min(np.abs(temps-t_fin)))
# Si les instants donnes sont exactement a mi-temsp entre deux instants de releve :
#    -> on choisi la plus petite valeur de debut et la plus grande valeur de fin
if (len(it_deb)==2) : it_deb = it_deb[0]
if (len(it_fin)==2) : it_fin = it_fin[-1]
# Si on a donne le meme instant de debut et de fin :
#   -> on demande les données à un seul instant.
#   --> /!\ Pour les vitesses, il faut prendre +2 pour pouvoir deriver numériquement
if (it_fin==it_deb) : 
    it_fin +=2
    # Si it_fin+2 > temps.shape[0]
    it_fin, it_deb = it_fin-(it_fin-temps.shape[0]), it_deb-(it_fin-temps.shape[0]) 
it_deb, it_fin = int(it_deb), int(it_fin)
print(len(temps))
print("# Debut des releves : "+str(temps[0])+' it :'+str(it_deb))
print("# Fin des releves : "+str(temps[-1])+'  it :'+str(it_fin))

########################################################################
##############################################################################
# Lecture de la géométrie du domaine
# Lx = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"uniform_domain_size_i")
# Ly = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"uniform_domain_size_j")
# Lz = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"uniform_domain_size_k")
# Ox = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"origin_i")
# Oy = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"origin_j")
# Oz = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"origin_k")
# Nx = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nbelem_i")     )
# Ny = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nbelem_j")     )
# Nz = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nbelem_k")     )
# Nt = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nb_pas_dt_max"))
# dx, dy, dz = Lx/Nx, Ly/Ny, Lz/Nz
# ex, ey, ez = Lx/float(Nx), Ly/float(Ny), Lz/float(Nz)

# Valeurs de references (a lire dans .data) : diametre bulle, surface bulle sph., longueur domaine
db = 1e-03 # diamètre bulle
S = 4*np.pi*(db/2.)**2 # Surface bulle sphérique
L = 3.17e-02 # Longueur domaine

####################################################################
# Ce qu'on va faire
do_the_chi = True
##############################################################################
color = {"x":"r","y":"g","z":"b"}
marker = {"x":"x","y":"+","z":"."}
linestyle = {"x":"solid","y":"dashed","z":"dotted"}
##############################################################################
####################################################################

"""
Etape 1 : Lecture et preparation des donnes des sondes
"""
#### Vel
Sx_ux_r, Sx_uy_r, Sx_uz_r = Donnees(), Donnees(), Donnees()
Sy_ux_r, Sy_uy_r, Sy_uz_r = Donnees(), Donnees(), Donnees()
Sz_ux_r, Sz_uy_r, Sz_uz_r = Donnees(), Donnees(), Donnees()
U_r = [[Grandeur(Sx_ux_r), Grandeur(Sx_uy_r), Grandeur(Sx_uz_r)],
       [Grandeur(Sy_ux_r), Grandeur(Sy_uy_r), Grandeur(Sy_uz_r)],
       [Grandeur(Sz_ux_r), Grandeur(Sz_uy_r), Grandeur(Sz_uz_r)]]
U = [U_r]
DIRECTION = [[r'$u_x$', r'$u_y$',r'$u_z$'],
             [r'$u_x$', r'$u_y$',r'$u_z$'],
             [r'$u_x$', r'$u_y$',r'$u_z$']]
for i,direction in enumerate(direction_sondes):
    for j,grandeur_et_direction in enumerate(direction_forces):
        velocity = grandeur_et_direction[-1]
        print("velocity",velocity)
        for u in U : 
            u[i][j].repoSondes = chemin_son
            u[i][j].it_deb = it_deb
            u[i][j].it_fin = it_fin
            u[i][j].formule = DIRECTION[i][j]
            u[i][j].liste_nomSondes = []
            u[i][j].nom_grandeur = 'force selon {0}'.format(velocity[-1]) 
            u[i][j].nom_sonde = grandeur_et_direction 
            u[i][j].nom_fichier = 'spectral_point2'
            u[i][j].direction_sonde = 'S_{0}'.format(direction)
            u[i][j].direction_vitesse = velocity
            u[i][j].identification = 'OO'
            # On a un signal tres tres bruite en temps, on invente donc la smoothing window
            u[i][j].smoothing_window = np.ones(N_smooth)
            u[i][j].get_the_n();
            u[i][j].initialize_datas()
            u[i][j].get_crude_from_probe()
            
        """ Cette partie va etre refaite, apres la relocalisation """
        for u in U:
            u[i][j].add_T=True
            u[i][j].add_S=True
            u[i][j].get_all_stats()
            u[i][j].get_all_stats_time()
            u[i][j].get_all_stats_space()
        U_r[i][j].get_regular_fluctu()
        U_r[i][j].get_regular_fluctu_time()
        U_r[i][j].get_regular_fluctu_space()
        for u in U:
            compute_correlations(u[i][j])
            # a coder : compute_correlations(u[i][j],t_deb,t_fin)
            
temps_corr = temps[0:int((len(temps)-1)/2)]-temps[0]
# Sorties - Temporel ###################################################

fig_Val, (ax_Val) = plt.subplots(1,1,figsize=(10,10))
fig_Val_smooth_mean, (ax_Val_smooth_mean) = plt.subplots(1,1,figsize=(10,10))
fig_CorT, (ax_CorT) = plt.subplots(1,1,figsize=(10,10))
fig_CorT2, (ax_CorT2) = plt.subplots(1,1,figsize=(10,10))
fig_CorT2_diag, (ax_CorT2_diag) = plt.subplots(1,1,figsize=(10,10))

fig_Val.suptitle("Valeurs sondes - Temporel",fontsize=24)
fig_Val_smooth_mean.suptitle("Valeurs sondes - Temporel - avg {0} elements".format(len(U[0][i][j].smoothing_window)),fontsize=24)
fig_CorT.suptitle(r"Autocorrelation $\frac{\langle %s_i %s_i \rangle}{\langle %s_i ^2 \rangle}$ sondes - Temporel"%(quelle_grandeur.lower(),quelle_grandeur.lower(),quelle_grandeur.lower()),fontsize=24)
fig_CorT2.suptitle(r"Autocorrelation $\frac{\langle %s_i %s_i \rangle}{\langle %s_i ^2 \rangle}$ sondes - Temporel - normal"%(quelle_grandeur.lower(),quelle_grandeur.lower(),quelle_grandeur.lower()),fontsize=24)
fig_CorT2_diag.suptitle(r"Autocorrelation $\frac{\langle %s_i %s_i \rangle}{\langle %s_i ^2 \rangle}$ sondes - Temporel - normal"%(quelle_grandeur.lower(),quelle_grandeur.lower(),quelle_grandeur.lower()),fontsize=24)

liste_label_sonde, liste_label_sonde_diag = [], []
liste_moyenne_sondes,liste_variance_sondes, liste_skewness_sondes, liste_kurtosis_sondes = [], [], [], []
liste_moyenne_sondes_norm,liste_variance_sondes_norm, liste_skewness_sondes_norm, liste_kurtosis_sondes_norm = [], [], [], []

for i,direction in enumerate(direction_sondes):
    for j,velocity in enumerate(direction_forces): 
        
        U_corr_T_norm = U[0][i][j].correlation_temporelle[int(len(U[0][i][j].correlation_temporelle)/2):-1] / U[0][i][j].correlation_temporelle[int(len(U[0][i][j].correlation_temporelle)/2)] 
        U_corr_T2_norm = U[0][i][j].correlation_temporelle2[int(len(U[0][i][j].correlation_temporelle2)/2):-1] / U[0][i][j].correlation_temporelle2[int(len(U[0][i][j].correlation_temporelle2)/2)] 
        
        ax_Val.plot(temps[1:],
                    U[0][i][j].fluctuations[2,:,int(U[0][i][j].npoints/2)],
                    linewidth=0.8,
                    color=color[velocity[-1].lower()],
                    linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
        ax_Val_smooth_mean.plot(temps[1:],
                                U[0][i][j].fluctuations_T2[2,:,int(U[0][i][j].npoints/2)],
                                linewidth=0.8,
                                color=color[velocity[-1].lower()],
                                linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))

        ax_CorT.plot(temps_corr[0:len(U_corr_T_norm)],
                     U_corr_T_norm,
                     linewidth=0.8,
                     color=color[velocity[-1].lower()],
                     linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
        ax_CorT2.plot(temps_corr[0:len(U_corr_T2_norm)],
                      U_corr_T2_norm,
                      linewidth=0.8,
                      color=color[velocity[-1].lower()],
                      linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
        # Uniquement les termes diagonaux
        if direction[0].lower()==velocity[-1].lower():
            # Le terme diagonal croise bien l'axe des abscisses
            if len(np.where(U_corr_T2_norm<=0)[0]) != 0:
                ax_CorT2_diag.plot(temps_corr[0:len(U_corr_T2_norm)],
                                    U_corr_T2_norm,
                                    linewidth=1.2,
                                    color=color[velocity[-1].lower()],
                                    linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
                ax_CorT2_diag.scatter(temps_corr[np.where(U_corr_T2_norm<=0)[0][0]],
                                    U_corr_T2_norm[np.where(U_corr_T2_norm<=0)[0][0]],
                                    color=color[velocity[-1].lower()])
                liste_label_sonde_diag.append(direction + U[0][i][j].nom_sonde[-2:] + r"; $\langle f_i ^2 \rangle =$ " + str(np.round(U[0][i][j].correlation_temporelle[int(len(U[0][i][j].correlation_temporelle)/2)],2)))
                                    
        liste_label_sonde.append(direction + U[0][i][j].nom_sonde[-2:] + r"; $\langle f_i ^2 \rangle =$ " + str(np.round(U[0][i][j].correlation_temporelle[int(len(U[0][i][j].correlation_temporelle)/2)],2)))
        # Liste des statistiques "par instant"
        liste_moyenne_sondes.append(U[0][i][j].moyenne_T)
        liste_variance_sondes.append(U[0][i][j].variance_T)
        liste_skewness_sondes.append(U[0][i][j].skewness_T)
        liste_kurtosis_sondes.append(U[0][i][j].kurtosis_T)
        # Liste des statistiques "conventionnelles"
        liste_moyenne_sondes_norm.append(U[0][i][j].moyenne)
        liste_variance_sondes_norm.append(U[0][i][j].variance)
        liste_skewness_sondes_norm.append(U[0][i][j].skewness)
        liste_kurtosis_sondes_norm.append(U[0][i][j].kurtosis)
        
ax_Val.legend(liste_label_sonde)    
ax_Val_smooth_mean.legend(liste_label_sonde)    
ax_CorT.legend(liste_label_sonde)    
ax_CorT2.legend(liste_label_sonde)    
ax_CorT2_diag.legend(liste_label_sonde_diag)    

ax_Val.set_xlabel('temps (s)',fontsize=24)
ax_Val_smooth_mean.set_xlabel('temps (s)',fontsize=24)
ax_CorT.set_xlabel('temps (s)',fontsize=24)
ax_CorT2.set_xlabel('temps (s)',fontsize=24)
ax_CorT2_diag.set_xlabel('temps (s)',fontsize=24)

ax_Val.set_ylabel(            r"$ %s_{ph}$"%(quelle_grandeur),fontsize=24)        
ax_Val_smooth_mean.set_ylabel(r"$ %s_{ph}$"%(quelle_grandeur),fontsize=24)        
ax_CorT.set_ylabel(           r"$\langle %s_{ph}; %s_{ph} \rangle$"%(quelle_grandeur,quelle_grandeur),fontsize=24)        
ax_CorT2.set_ylabel(          r"$\langle %s_{ph}; %s_{ph} \rangle$"%(quelle_grandeur,quelle_grandeur),fontsize=24)        
ax_CorT2_diag.set_ylabel(     r"$\langle %s_{ph}; %s_{ph} \rangle$"%(quelle_grandeur,quelle_grandeur),fontsize=24)        

ax_Val.grid(True,which='both')        
ax_Val_smooth_mean.grid(True,which='both')        
ax_CorT.grid(True,which='both')        
ax_CorT2.grid(True,which='both')        
ax_CorT2_diag.grid(True,which='both')        

fig_Val.savefig(quelle_grandeur+"_"+"Force_T.png")     
fig_Val_smooth_mean.savefig(quelle_grandeur+"_"+"Force_T_mean.png")     
fig_CorT.savefig(quelle_grandeur+"_"+"Correlation_T.png")     
fig_CorT2.savefig(quelle_grandeur+"_"+"Correlation_T2.png")     
fig_CorT2_diag.savefig(quelle_grandeur+"_"+"Correlation_T2_diagonaux.png")     

print("##############################################")
print("# Statistiques des sondes relevant la force spectrale")
# print("# Statistiques \"par point de l'espace\"")
# print("# Liste des sondes : "+str(liste_label_sonde))
# print("# Liste des moyennes : "+str(liste_moyenne_sondes))
# print("# Liste des variances : "+str(liste_variance_sondes))
# print("# Liste des skewness (moment d'o 3) : "+str(liste_skewness_sondes))
# print("# Liste des kurtosis (moment d'o 4) : "+str(liste_kurtosis_sondes))
print("####")      
print("# Statistiques \"conventionnelles\" des sondes")
print(" sondes : "+str(liste_label_sonde))
print("\"moyennes\" : "+str(liste_moyenne_sondes_norm)+",")
print("\"variances\" : "+str(liste_variance_sondes_norm)+",")
print("\"skewness\" : "+str(liste_skewness_sondes_norm)+",")
print("\"kurtosis\" : "+str(liste_kurtosis_sondes_norm))
print("}")
print("####")
print("# Statistiques \"conventionnelles\" des sondes  - arrondis")
lcm = np.floor(np.log10(np.abs(np.array(liste_moyenne_sondes_norm))))
lcv = np.floor(np.log10(np.abs(np.array(liste_variance_sondes_norm))))
lcs = np.floor(np.log10(np.abs(np.array(liste_skewness_sondes_norm))))
lck = np.floor(np.log10(np.abs(np.array(liste_kurtosis_sondes_norm))))
liste_moyenne_sondes_norm = np.round(liste_moyenne_sondes_norm*10**(-lcm),3)
liste_variance_sondes_norm = np.round(liste_variance_sondes_norm*10**(-lcv),3)
liste_skewness_sondes_norm = np.round(liste_skewness_sondes_norm*10**(-lcs),3)
liste_kurtosis_sondes_norm = np.round(liste_kurtosis_sondes_norm*10**(-lck),3)
print("# Liste des sondes : "+str(liste_label_sonde)                      )
print("# (ca) Liste des expo moy : "+str(lcm)                 )
print("# (ca) Liste des moyennes : "+str(liste_moyenne_sondes_norm)                 )
print("# (ca) Liste des expo var : "+str(lcv)               )
print("# (ca) Liste des variances : "+str(liste_variance_sondes_norm)               )
print("# (ca) Liste des expo ske (moment d'o 3) : "+str(lcs) )
print("# (ca) Liste des skewness (moment d'o 3) : "+str(liste_skewness_sondes_norm) )
print("# (ca) Liste des expo kur (moment d'o 4) : "+str(lck) )
print("# (ca) Liste des kurtosis (moment d'o 4) : "+str(liste_kurtosis_sondes_norm) )
print("##############################################")

                
########################################################################                
# Sorties - Spatial ####################################################

fig_Val, (ax_Val) = plt.subplots(1,1,figsize=(10,10))
fig_CorS, (ax_CorS) = plt.subplots(1,1,figsize=(10,10))
fig_CorS2, (ax_CorS2) = plt.subplots(1,1,figsize=(10,10))
fig_CorS2_diag, (ax_CorS2_diag) = plt.subplots(1,1,figsize=(10,10))
fig_Val.suptitle("Valeurs sondes - Spatial",fontsize=24)
fig_CorS.suptitle(r"Autocorrelation sondes $\frac{\langle %s_i %s_i \rangle}{\langle %s_i ^2 \rangle}$ - Spatial"%(quelle_grandeur.lower(),quelle_grandeur.lower(),quelle_grandeur.lower()),fontsize=24)
fig_CorS2.suptitle(r"Autocorrelation sondes $\frac{\langle %s_i %s_i \rangle}{\langle %s_i ^2 \rangle}$ - Spatial - normal"%(quelle_grandeur.lower(),quelle_grandeur.lower(),quelle_grandeur.lower()),fontsize=24)
fig_CorS2_diag.suptitle(r"Autocorrelation sondes $\frac{\langle %s_i %s_i \rangle}{\langle %s_i ^2 \rangle}$ - Spatial - normal"%(quelle_grandeur.lower(),quelle_grandeur.lower(),quelle_grandeur.lower()),fontsize=24)

liste_label_sonde, liste_label_sonde_diag = [], []
liste_moyenne_sondes, liste_variance_sondes, liste_skewness_sondes, liste_kurtosis_sondes = [], [], [], []
liste_moyenne_sondes_norm, liste_variance_sondes_norm, liste_skewness_sondes_norm, liste_kurtosis_sondes_norm = [], [], [], []

espace = U[0][i][j].coord[direction_sondes[i][0].swapcase()]
espace_corr = espace[2:int((len(espace)-1)/2)+2]; n=len(espace_corr)
espace_corr = espace_corr - espace_corr[2]
# espace_corr = espace_corr[direction[0].swapcase()][-1]-espace_corr[direction[0].swapcase()]
print(espace.shape)
print(U[0][i][j].correlation_spatiale.shape)
print(U[0][i][j].correlation_spatiale2.shape)
print(U[0][i][j].fluctuations_S.shape)
print(U[0][i][j].fluctuations_S2.shape)

for i,direction in enumerate(direction_sondes):
    for j,velocity in enumerate(direction_forces): 
        
        U_corr_S_norm = (U[0][i][j].correlation_spatiale[int(len(U[0][i][j].correlation_spatiale)/2):-1] / U[0][i][j].correlation_spatiale[int(len(U[0][i][j].correlation_spatiale)/2)])
        U_corr_S2_norm = (U[0][i][j].correlation_spatiale2[int(len(U[0][i][j].correlation_spatiale2)/2):-1] / U[0][i][j].correlation_spatiale2[int(len(U[0][i][j].correlation_spatiale2)/2)])
        
        ax_Val.plot(espace,
                    U[0][i][j].fluctuations[2,int(U[0][i][j].ntemps/2),:],
                    linewidth=0.8,
                    color=color[velocity[-1].lower()],
                    linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
        # ax_CorS.scatter(espace[1:]-espace[len(espace)-1]/4.,U[0][i][j].correlation_spatiale,linewidth=0.8,color=color[j],marker=marker[i])#,label=str(label_sonde))
        ax_CorS.plot(espace_corr,
                     U_corr_S_norm,
                     linewidth=0.8,
                     color=color[velocity[-1].lower()],
                     linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
        ax_CorS2.plot(espace_corr,
                      U_corr_S2_norm,
                      linewidth=0.8,
                      color=color[velocity[-1].lower()],
                      linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
        # Uniquement les temres diagonaux
        if velocity[-1].lower()==direction[0].lower():
            # Le terme diagonal croise bien l'axe des abscisses
            if len(np.where(U_corr_S2_norm<=0)[0]) != 0:
                ax_CorS2_diag.plot(espace_corr,
                                   U_corr_S2_norm,
                                   linewidth=1.2,
                                   color=color[velocity[-1].lower()],
                                   linestyle=linestyle[direction[0].lower()])#,label=str(label_sonde))
                ax_CorS2_diag.scatter(espace_corr[2+np.where(U_corr_S2_norm<=0)[0][0]],
                                      U_corr_S2_norm[np.where(U_corr_S2_norm<=0)[0][0]],
                                      color=color[velocity[-1].lower()])
                liste_label_sonde_diag.append(direction + U[0][i][j].nom_sonde[-2:] + r"; $\langle f_i ^2 \rangle =$ " + str(np.round(U[0][i][j].correlation_spatiale2[int(len(U[0][i][j].correlation_spatiale2)/2)],3)))

        liste_label_sonde.append(direction + U[0][i][j].nom_sonde[-2:] + r"; $\langle f_i ^2 \rangle =$ " + str(np.round(U[0][i][j].correlation_spatiale2[int(len(U[0][i][j].correlation_spatiale2)/2)],3)))
        # Liste des statistiques "par instant"
        liste_moyenne_sondes.append(U[0][i][j].moyenne_S)
        liste_variance_sondes.append(U[0][i][j].variance_S)
        liste_skewness_sondes.append(U[0][i][j].skewness_S)
        liste_kurtosis_sondes.append(U[0][i][j].kurtosis_S)
        # Liste des statistiques "conventionnelles "
        liste_moyenne_sondes_norm.append(U[0][i][j].moyenne)
        liste_variance_sondes_norm.append(U[0][i][j].variance)
        liste_skewness_sondes_norm.append(U[0][i][j].skewness)
        liste_kurtosis_sondes_norm.append(U[0][i][j].kurtosis)
        
ax_Val.legend(liste_label_sonde)    
ax_CorS.legend(liste_label_sonde)    
ax_CorS2.legend(liste_label_sonde)    
ax_CorS2_diag.legend(liste_label_sonde_diag)    
ax_Val.set_xlabel('x, y ou z (m)',fontsize=24)
ax_CorS.set_xlabel('x,y ou z (m)',fontsize=24)
ax_CorS2.set_xlabel('x,y ou z (m)',fontsize=24)
ax_CorS2_diag.set_xlabel('x,y ou z (m)',fontsize=24)
ax_Val.set_ylabel(       r"$ %s_{ph}$"%(quelle_grandeur),fontsize=24)        
ax_CorS.set_ylabel(      r"$\langle %s_{ph}; %s_{ph} \rangle$"%(quelle_grandeur,quelle_grandeur),fontsize=24)        
ax_CorS2.set_ylabel(     r"$\langle %s_{ph}; %s_{ph} \rangle$"%(quelle_grandeur,quelle_grandeur),fontsize=24)        
ax_CorS2_diag.set_ylabel(r"$\langle %s_{ph}; %s_{ph} \rangle$"%(quelle_grandeur,quelle_grandeur),fontsize=24)        
ax_Val.grid(True,which='both')        
ax_CorS.grid(True,which='both')        
ax_CorS2.grid(True,which='both')        
ax_CorS2_diag.grid(True,which='both')        
fig_Val.savefig(quelle_grandeur+"_"+"Force_S.png")     
fig_CorS.savefig(quelle_grandeur+"_"+"Correlation_S.png")     
fig_CorS2.savefig(quelle_grandeur+"_"+"Correlation_S2.png")     
fig_CorS2_diag.savefig(quelle_grandeur+"_"+"Correlation_S2_diagonaux.png")     

print("##############################################")
# print("# Statistiques des sondes \"par instant\"")
# print("# Liste des sondes : "+str(liste_label_sonde))
# print("# Liste des moyennes : "+str(liste_moyenne_sondes))
# print("# Liste des variances : "+str(liste_variance_sondes))
# print("# Liste des skewness (moment d'o 3) : "+str(liste_skewness_sondes))
# print("# Liste des kurtosis (moment d'o 4) : "+str(liste_kurtosis_sondes))
print("####")
print("# Statistiques conventionnelles des sondes (c)")
print("sondes : "+str(liste_label_sonde)+",")
print("{")
print("\"moyennes\" : "+str(liste_moyenne_sondes_norm)+",")
print("\"variances\" : "+str(liste_variance_sondes_norm)+",")
print("\"skewness\" : "+str(liste_skewness_sondes_norm)+",")
print("\"kurtosis\" : "+str(liste_kurtosis_sondes_norm))
print("}")
print("####")
print("# Statistiques des sondes \"conventionnelles\" - arrondis")
print(np.log10(np.abs(np.array(liste_moyenne_sondes))))
lcm = np.floor(np.log10(np.abs(np.array(liste_moyenne_sondes_norm))))
lcv = np.floor(np.log10(np.abs(np.array(liste_variance_sondes_norm))))
lcs = np.floor(np.log10(np.abs(np.array(liste_skewness_sondes_norm))))
lck = np.floor(np.log10(np.abs(np.array(liste_kurtosis_sondes_norm))))
liste_moyenne_sondes_norm = np.round(liste_moyenne_sondes_norm*10**(-lcm),3)
liste_variance_sondes_norm = np.round(liste_moyenne_sondes_norm*10**(-lcv),3)
liste_skewness_sondes_norm = np.round(liste_moyenne_sondes_norm*10**(-lcs),3)
liste_kurtosis_sondes_norm = np.round(liste_moyenne_sondes_norm*10**(-lck),3)
print("# Liste des sondes : "+str(liste_label_sonde)                      )
print("# (ca) Liste des expo moy : "+str(lcm)                 )
print("# (ca) Liste des moyennes : "+str(liste_moyenne_sondes_norm)                 )
print("# (ca) Liste des expo var : "+str(lcv)               )
print("# (ca) Liste des variances : "+str(liste_variance_sondes_norm)               )
print("# (ca) Liste des expo ske (moment d'o 3) : "+str(lcs) )
print("# (ca) Liste des skewness (moment d'o 3) : "+str(liste_skewness_sondes_norm) )
print("# (ca) Liste des expo kur (moment d'o 4) : "+str(lck) )
print("# (ca) Liste des kurtosis (moment d'o 4) : "+str(liste_kurtosis_sondes_norm) )
print("##############################################")


########################################################################                
# Obtention échelles #######################################
# les prints sont pensés pour que le fichier de sortie soit recopiable 
# en python de sorte à definir des dictionnaires, facilement
# integartion : simps chez scipy mieux que mon truc
for i,direction in enumerate(direction_sondes):
    for j,velocity in enumerate(direction_forces): 
        if velocity[-1].lower()==direction[0].lower():  
            U_corr_S2_norm = (U[0][i][j].correlation_spatiale2[int(len(U[0][i][j].correlation_spatiale2)/2):-1] / U[0][i][j].correlation_spatiale2[int(len(U[0][i][j].correlation_spatiale2)/2)])
            U_corr_T2_norm = U[0][i][j].correlation_temporelle2[int(len(U[0][i][j].correlation_temporelle2)/2):-1] / U[0][i][j].correlation_temporelle2[int(len(U[0][i][j].correlation_temporelle2)/2)] 
            
            if len(np.where(U_corr_T2_norm<=0)[0]) != 0:
                integ_T= simps(U_corr_T2_norm[0:np.where(U_corr_T2_norm<=0)[0][0]],temps_corr[0:np.where(U_corr_T2_norm<=0)[0][0]])
            if len(np.where(U_corr_S2_norm<=0)[0]) != 0:
                integ_S= simps(U_corr_S2_norm[0:np.where(U_corr_S2_norm<=0)[0][0]],espace_corr[2:2+np.where(U_corr_S2_norm<=0)[0][0]])
            
            ## f(x+dx) -2f(x) + f(x-dx) = dx**2 f''(x)
            # ==> pour f paire et pour x=0 : 2 (f(dx)-f(0)) = dx**2 f''(0)
            # ==> L_taylor = (-0.5 f''(0)) ** (-1/2)
            Taylor_S = (-(2./(espace_corr[3]-espace_corr[2])**2) * (U_corr_S2_norm[1]-U_corr_S2_norm[0]))**(-1./2)
            Taylor_T = (-(2./(temps_corr[1]-temps_corr[0])  **2) * (U_corr_T2_norm[1]-U_corr_T2_norm[0]))**(-1./2)
            
            print("####################################################")
            print("### Grandeur caractéristiques")
            print("\""+velocity[-1].lower()+"\" : {")
            #
            print("\"Annulation de l'autocorrelation normée\" : ")
            if len(np.where(U_corr_T2_norm<=0)[0]) != 0:
                print("{\"T\" : "+str(temps_corr[np.where(U_corr_T2_norm<=0)[0][0]])+",")
            else :
                print("{\"T\" : \"L'autocorrelation normee ne s'annule pas...\",")
            if len(np.where(U_corr_S2_norm<=0)[0]) != 0:
                print("\"L\" : "+str(espace_corr[2+np.where(U_corr_S2_norm<=0)[0][0]])+"},")
            else :
                print("\"L\" : \"L'autocorrelation normee ne s'annule pas...\"},")
            #
            print("\"Integrale entre 0 et annulation\" : ")
            if len(np.where(U_corr_T2_norm<=0)[0]) != 0:
                print("{\"T\" : "+str(integ_T)+",")
            else :
                print("{\"T\" : \"L'autocorrelation normee ne s'annule pas...\",")    
            if len(np.where(U_corr_S2_norm<=0)[0]) != 0:
                print("\"L\" : "+str(integ_S)+"},")
            else :
                print("\"L\" : \"L'autocorrelation normee ne s'annule pas...\"},")
            #
            print("\"Micro-échelle de Taylor (definie pour L, transposée pour T)\" : ")
            print("{\"L\" : "+str(Taylor_S)+",")
            print("\"T\" : "+str(Taylor_T)+"},")
            print("},")
            
# 
