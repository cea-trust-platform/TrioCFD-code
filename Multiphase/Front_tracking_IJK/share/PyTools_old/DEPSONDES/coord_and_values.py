# -*- coding: utf8

import os
import numpy as np
import subprocess
import scipy
from scipy import signal
import pandas
import glob
import sys

# import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use("Agg")
# import matplotlib.animation as animation

from depSondes_classes_et_fonctions import *
import DNSTools3 as dtool

"""
    Objectif : tracer les positions de toutes les sondes du dossier SON sondes
               tracer les valeurs de toutes les grandeurs des sondes du dossier SON
"""

##############################################################################
# mode d'emploi : python coord_and_values.py nom_data chemin_son liste_grandeurs
#                                              1         2            3
# /!\ 1 : ne doit PAS terminer par .data
# /!\ 3 : [code_grandeur1,code_grandeur2, ...] -> voir la variable "rosette_noms"
#         dans .../GENSONDES/genSondes_GR.py
##############################################################################
liste_grandeurs = sys.argv[3].lstrip("[").rstrip("]").split(",")
print(liste_grandeurs)
chemin_son = sys.argv[2]
nom_data = sys.argv[1]
##############################################################################
# Lecture de la géométrie du domaine
Lx = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"uniform_domain_size_i")
Ly = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"uniform_domain_size_j")
Lz = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"uniform_domain_size_k")
Ox = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"origin_i")
Oy = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"origin_j")
Oz = dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"origin_k")
Nx = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nbelem_i")     )
Ny = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nbelem_j")     )
Nz = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nbelem_k")     )
Nt = int(dtool.getParam(glob.glob(chemin_son+"../"+nom_data+".data")[0],"nb_pas_dt_max"))
dx, dy, dz = Lx/Nx, Ly/Ny, Lz/Nz
ex, ey, ez = Lx/float(Nx), Ly/float(Ny), Lz/float(Nz)
##############################################################################
color = ["r","g","b"]
marker = ["x","+","."]
##############################################################################

##############################################################################
#   COORDONNEES
##############################################################################

liste_label_sonde = []
fig_XY, (ax_XY) = plt.subplots(1,1,figsize=(10,10))
fig_XZ, (ax_XZ) = plt.subplots(1,1,figsize=(10,10))
fig_YZ, (ax_YZ) = plt.subplots(1,1,figsize=(10,10))

fig_XY.suptitle("Positions sondes - Plan XY",fontsize=24)
fig_XZ.suptitle("Positions sondes - Plan XZ",fontsize=24)
fig_YZ.suptitle("Positions sondes - Plan YZ",fontsize=24)

ax_XY.set_xlim(Ox-ex,Ox+Lx+ex)
ax_XY.set_ylim(Oy-ey,Oy+Ly+ey)
ax_XZ.set_xlim(Ox-ex,Ox+Lx+ex)
ax_XZ.set_ylim(Oz-ez,Oz+Lz+ez)
ax_YZ.set_xlim(Oy-ey,Oy+Ly+ey)
ax_YZ.set_ylim(Oz-ez,Oz+Lz+ez)

for iS, direction_S in enumerate(["_S_X","_S_Y","_S_Z"]):
    for iV, direction_V in enumerate(["X","Y","Z"]):
        liste_nom_sonde = glob.glob(chemin_son+nom_data+direction_S+"*11*"+direction_V+".son")
        for sonde in liste_nom_sonde:
            label_sonde = sonde.replace(chemin_son+nom_data,"").lstrip("_")
            label_sonde = label_sonde.rstrip(".son")
            liste_label_sonde.append(label_sonde)
            # print(sonde)
            print(label_sonde)
            coord = coordSonde(sonde)
            # print(coord)
            ax_XY.scatter(coord["x"],coord["y"],color=color[iV],marker=marker[iS])#,label=str(label_sonde))
            ax_XZ.scatter(coord["x"],coord["z"],color=color[iV],marker=marker[iS])#,label=str(label_sonde))
            ax_YZ.scatter(coord["y"],coord["z"],color=color[iV],marker=marker[iS])#,label=str(label_sonde))
            
ax_XY.legend(liste_label_sonde)    
ax_XZ.legend(liste_label_sonde)    
ax_YZ.legend(liste_label_sonde)    

ax_XY.set_xlabel('X',fontsize=24)
ax_XZ.set_xlabel('X',fontsize=24)
ax_YZ.set_xlabel('Y',fontsize=24)

ax_XY.set_ylabel('Y',fontsize=24)        
ax_XZ.set_ylabel('Z',fontsize=24)        
ax_YZ.set_ylabel('Z',fontsize=24) 

ax_XY.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_XZ.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_YZ.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)

ax_XY.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_XZ.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_YZ.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
    
ax_XY.set_yticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_XZ.set_yticks(np.arange(Oz,Oz+Lz,dz),minor=True)
ax_YZ.set_yticks(np.arange(Oz,Oz+Lz,dz),minor=True)  

ax_XY.set_yticks(np.arange(Oy,Oy+Ly,5*dy))
ax_XZ.set_yticks(np.arange(Oz,Oz+Lz,5*dz))
ax_YZ.set_yticks(np.arange(Oz,Oz+Lz,5*dz))       

ax_XY.grid(True,which='both')        
ax_XZ.grid(True,which='both')        
ax_YZ.grid(True,which='both')        

ax_XY.grid(True,which='minor', alpha=0.2, linestyle='--')        
ax_XZ.grid(True,which='minor', alpha=0.2, linestyle='--')        
ax_YZ.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 

fig_XY.savefig("coord_XY.png")           
fig_XZ.savefig("coord_XZ.png")           
fig_YZ.savefig("coord_YZ.png")           
            
##############################################################################
#   VALEURS
##############################################################################

fig_X_init, (ax_X_init) = plt.subplots(1,1,figsize=(10,10))
fig_Y_init, (ax_Y_init) = plt.subplots(1,1,figsize=(10,10))
fig_Z_init, (ax_Z_init) = plt.subplots(1,1,figsize=(10,10))

fig_X_med, (ax_X_med) = plt.subplots(1,1,figsize=(10,10))
fig_Y_med, (ax_Y_med) = plt.subplots(1,1,figsize=(10,10))
fig_Z_med, (ax_Z_med) = plt.subplots(1,1,figsize=(10,10))

fig_X_last, (ax_X_last) = plt.subplots(1,1,figsize=(10,10))
fig_Y_last, (ax_Y_last) = plt.subplots(1,1,figsize=(10,10))
fig_Z_last, (ax_Z_last) = plt.subplots(1,1,figsize=(10,10))

fig_X_init.suptitle("Valeurs sondes - Plan XY",fontsize=24)
fig_Y_init.suptitle("Valeurs sondes - Plan XZ",fontsize=24)
fig_Z_init.suptitle("Valeurs sondes - Plan YZ",fontsize=24)


fig_X_med.suptitle("Valeurs sondes - Plan XY",fontsize=24)
fig_Y_med.suptitle("Valeurs sondes - Plan XZ",fontsize=24)
fig_Z_med.suptitle("Valeurs sondes - Plan YZ",fontsize=24)


fig_X_last.suptitle("Valeurs sondes - Plan XY",fontsize=24)
fig_Y_last.suptitle("Valeurs sondes - Plan XZ",fontsize=24)
fig_Z_last.suptitle("Valeurs sondes - Plan YZ",fontsize=24)

ax_X_init.set_xlim(Ox-ex,Ox+Lx+ex)
ax_X_init.set_xlim(Ox-ex,Ox+Lx+ex)
ax_X_init.set_xlim(Oy-ey,Oy+Ly+ey)

ax_Y_init.set_xlim(Ox-ex,Ox+Lx+ex)
ax_Y_init.set_xlim(Ox-ex,Ox+Lx+ex)
ax_Y_init.set_xlim(Oy-ey,Oy+Ly+ey)

ax_Z_init.set_xlim(Ox-ex,Ox+Lx+ex)
ax_Z_init.set_xlim(Ox-ex,Ox+Lx+ex)
ax_Z_init.set_xlim(Oy-ey,Oy+Ly+ey)

for iS, direction_S in enumerate(["_S_X","_S_Y","_S_Z"]):
    for iV, direction_V in enumerate(["X","Y","Z"]):
        for iG, grandeur in enumerate(liste_grandeurs) :
            liste_nom_sonde = glob.glob(chemin_son+nom_data+direction_S+"*11*"+"_"+grandeur+direction_V+".son")
            for sonde in liste_nom_sonde:
                label_sonde = sonde.replace(chemin_son+nom_data,"").lstrip("_")
                label_sonde = label_sonde.rstrip(".son")
                liste_label_sonde.append(label_sonde)
                print(label_sonde)
                coord = coordSonde(sonde)
                
                valeur = loadSondeSeg(sonde)
                n = valeur.shape[0]; n_med = int(n/2)
                print(n)
                
                ax_X_init.scatter(coord["x"],valeur[0][1:],color=color[iV],
                                marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
                ax_Y_init.scatter(coord["y"],valeur[0][1:],color=color[iV],
                                marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
                ax_Z_init.scatter(coord["z"],valeur[0][1:],color=color[iV],
                             marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
            
                    
                ax_X_med.scatter(coord["x"],valeur[n_med][1:],color=color[iV],
                            marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
                ax_Y_med.scatter(coord["y"],valeur[n_med][1:],color=color[iV],
                            marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
                ax_Z_med.scatter(coord["z"],valeur[n_med][1:],color=color[iV],
                            marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
            
                                
                ax_X_last.scatter(coord["x"],valeur[-1][1:],color=color[iV],
                            marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
                ax_Y_last.scatter(coord["y"],valeur[-1][1:],color=color[iV],
                            marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
                ax_Z_last.scatter(coord["z"],valeur[-1][1:],color=color[iV],
                            marker=marker[iS],s=40*(iG+1)**2)#,label=str(label_sonde))
            
# Legend
ax_X_init.legend(liste_label_sonde)    
ax_Y_init.legend(liste_label_sonde)    
ax_Z_init.legend(liste_label_sonde)    

ax_X_med.legend(liste_label_sonde)    
ax_Y_med.legend(liste_label_sonde)    
ax_Z_med.legend(liste_label_sonde)    

ax_X_last.legend(liste_label_sonde)    
ax_Y_last.legend(liste_label_sonde)    
ax_Z_last.legend(liste_label_sonde)    

# xlabel
ax_X_init.set_xlabel('X',fontsize=24)
ax_Y_init.set_xlabel('Y',fontsize=24)
ax_Z_init.set_xlabel('Z',fontsize=24)
    
ax_X_med.set_xlabel('X',fontsize=24)
ax_Y_med.set_xlabel('Y',fontsize=24)
ax_Z_med.set_xlabel('Z',fontsize=24)
    
ax_X_last.set_xlabel('X',fontsize=24)
ax_Y_last.set_xlabel('Y',fontsize=24)
ax_Z_last.set_xlabel('Z',fontsize=24)

# ylabel
ax_X_init.set_ylabel(r"$F_{ph}$",fontsize=24)        
ax_Y_init.set_ylabel(r"$F_{ph}$",fontsize=24)        
ax_Z_init.set_ylabel(r"$F_{ph}$",fontsize=24)        
    
ax_X_med.set_ylabel(r"$F_{ph}$",fontsize=24)        
ax_Y_med.set_ylabel(r"$F_{ph}$",fontsize=24)        
ax_Z_med.set_ylabel(r"$F_{ph}$",fontsize=24)        
    
ax_X_last.set_ylabel(r"$F_{ph}$",fontsize=24)        
ax_Y_last.set_ylabel(r"$F_{ph}$",fontsize=24)        
ax_Z_last.set_ylabel(r"$F_{ph}$",fontsize=24)        

## grid
# For the minor
ax_X_init.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_X_init.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_X_init.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_X_med.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_X_med.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_X_med.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_X_last.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_X_last.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_X_last.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)

ax_Y_init.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Y_init.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Y_init.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_Y_med.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Y_med.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Y_med.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_Y_last.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Y_last.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Y_last.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)

ax_Z_init.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Z_init.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Z_init.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_Z_med.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Z_med.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Z_med.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)
ax_Z_last.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Z_last.set_xticks(np.arange(Ox,Ox+Lx,dx),minor=True)
ax_Z_last.set_xticks(np.arange(Oy,Oy+Ly,dy),minor=True)

# For the major
ax_X_init.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_X_init.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_X_init.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
ax_X_med.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_X_med.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_X_med.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
ax_X_last.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_X_last.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_X_last.set_xticks(np.arange(Oy,Oy+Ly,5*dy))

ax_Y_init.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Y_init.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Y_init.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
ax_Y_med.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Y_med.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Y_med.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
ax_Y_last.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Y_last.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Y_last.set_xticks(np.arange(Oy,Oy+Ly,5*dy))

ax_Z_init.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Z_init.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Z_init.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
ax_Z_med.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Z_med.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Z_med.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
ax_Z_last.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Z_last.set_xticks(np.arange(Ox,Ox+Lx,5*dx))
ax_Z_last.set_xticks(np.arange(Oy,Oy+Ly,5*dy))
    
# ax_XY.set_yticks(np.arange(Oy,Oy+Ly,dy),minor=True)
# ax_XZ.set_yticks(np.arange(Oz,Oz+Lz,dz),minor=True)
# ax_YZ.set_yticks(np.arange(Oz,Oz+Lz,dz),minor=True)  

# ax_XY.set_yticks(np.arange(Oy,Oy+Ly,5*dy))
# ax_XZ.set_yticks(np.arange(Oz,Oz+Lz,5*dz))
# ax_YZ.set_yticks(np.arange(Oz,Oz+Lz,5*dz))       

ax_X_init.grid(True,which='both')        
ax_X_init.grid(True,which='both')        
ax_X_init.grid(True,which='both')        
ax_X_med.grid(True,which='both')        
ax_X_med.grid(True,which='both')        
ax_X_med.grid(True,which='both')        
ax_X_last.grid(True,which='both')        
ax_X_last.grid(True,which='both')        
ax_X_last.grid(True,which='both')        

ax_X_init.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_X_init.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_X_init.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 
ax_X_med.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_X_med.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_X_med.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 
ax_X_last.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_X_last.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_X_last.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 

ax_Y_init.grid(True,which='both')        
ax_Y_init.grid(True,which='both')        
ax_Y_init.grid(True,which='both')        
ax_Y_med.grid(True,which='both')        
ax_Y_med.grid(True,which='both')        
ax_Y_med.grid(True,which='both')        
ax_Y_last.grid(True,which='both')        
ax_Y_last.grid(True,which='both')        
ax_Y_last.grid(True,which='both')        
   
ax_Y_init.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Y_init.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Y_init.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 
ax_Y_med.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Y_med.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Y_med.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 
ax_Y_last.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Y_last.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Y_last.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 

ax_Z_init.grid(True,which='both')        
ax_Z_init.grid(True,which='both')        
ax_Z_init.grid(True,which='both')        
ax_Z_med.grid(True,which='both')        
ax_Z_med.grid(True,which='both')        
ax_Z_med.grid(True,which='both')        
ax_Z_last.grid(True,which='both')        
ax_Z_last.grid(True,which='both')        
ax_Z_last.grid(True,which='both')        
   
ax_Z_init.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Z_init.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Z_init.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 
ax_Z_med.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Z_med.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Z_med.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 
ax_Z_last.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Z_last.grid(True,which='minor', alpha=0.2, linestyle='--')       
ax_Z_last.grid(True,which='minor', alpha=0.2, linestyle='--'    ) 


# Sauvegarde
fig_X_init.savefig("valeur_X_init.png")           
fig_Y_init.savefig("valeur_Y_init.png")           
fig_Z_init.savefig("valeur_Z_init.png")           
            
fig_X_med.savefig("valeur_X_med.png")           
fig_Y_med.savefig("valeur_Y_med.png")           
fig_Z_med.savefig("valeur_Z_med.png")           
     
fig_X_last.savefig("valeur_X_last.png")           
fig_Y_last.savefig("valeur_Y_last.png")           
fig_Z_last.savefig("valeur_Z_last.png")           
