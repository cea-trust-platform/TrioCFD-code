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

"""
    Objectif : tracer les positions des sondes
               tracer les valeurs des sondes
"""

##############################################################################
# mode d'emploi : python coord_and_velocity.py chemin_son
##############################################################################
chemin_son = sys.argv[1]
# ###############################################################
# Définition des dictionnaires
velocity = {"x":{"sx":"","sy":"","sz":""},
            "y":{"sx":"","sy":"","sz":""},
            "z":{"sx":"","sy":"","sz":""}
           }
coord_velocity = {"x":{"sx":0,"sy":0,"sz":0},
                  "y":{"sx":0,"sy":0,"sz":0},
                  "z":{"sx":0,"sy":0,"sz":0}
                 }           
valeur_velocity = {"x":{"sx":0,"sy":0,"sz":0},
                   "y":{"sx":0,"sy":0,"sz":0},
                   "z":{"sx":0,"sy":0,"sz":0}
                  }   
centered_velocity = {"x":{"sx":"","sy":"","sz":""},
                     "y":{"sx":"","sy":"","sz":""},
                     "z":{"sx":"","sy":"","sz":""}
                    }
coord_centered_velocity = {"x":{"sx":0,"sy":0,"sz":0},
                          "y":{"sx":0,"sy":0,"sz":0},
                          "z":{"sx":0,"sy":0,"sz":0}
                          }                    
valeur_centered_velocity = {"x":{"sx":0,"sy":0,"sz":0},
                            "y":{"sx":0,"sy":0,"sz":0},
                            "z":{"sx":0,"sy":0,"sz":0}
                            }                    
indicatrice = {"sx":"","sy":"","sz":""}
coord_indicatrice = {"sx":0,"sy":0,"sz":0}
valeur_indicatrice = {"sx":0,"sy":0,"sz":0}
# ##############################################################


# Chemins d'appel : grandeur[direction_grandeur][direction_sonde]
velocity["x"]["sx"] = chemin_son+"sondes_S_X_VELOCITY_X.son"
velocity["x"]["sy"] = chemin_son+"sondes_S_Y_VELOCITY_X.son"
velocity["x"]["sz"] = chemin_son+"sondes_S_Z_VELOCITY_X.son"

velocity["y"]["sx"] = chemin_son+"sondes_S_X_VELOCITY_Y.son"
velocity["y"]["sy"] = chemin_son+"sondes_S_Y_VELOCITY_Y.son"
velocity["y"]["sz"] = chemin_son+"sondes_S_Z_VELOCITY_Y.son"

velocity["z"]["sx"] = chemin_son+"sondes_S_X_VELOCITY_Z.son"
velocity["z"]["sy"] = chemin_son+"sondes_S_Y_VELOCITY_Z.son"
velocity["z"]["sz"] = chemin_son+"sondes_S_Z_VELOCITY_Z.son"

centered_velocity["x"]["sx"] = chemin_son+"sondes_S_X_CELL_VELOCITY_X.son"
centered_velocity["x"]["sy"] = chemin_son+"sondes_S_Y_CELL_VELOCITY_X.son"
centered_velocity["x"]["sz"] = chemin_son+"sondes_S_Z_CELL_VELOCITY_X.son"
                 
centered_velocity["y"]["sx"] = chemin_son+"sondes_S_X_CELL_VELOCITY_Y.son"
centered_velocity["y"]["sy"] = chemin_son+"sondes_S_Y_CELL_VELOCITY_Y.son"
centered_velocity["y"]["sz"] = chemin_son+"sondes_S_Z_CELL_VELOCITY_Y.son"
                 
centered_velocity["z"]["sx"] = chemin_son+"sondes_S_X_CELL_VELOCITY_Z.son"
centered_velocity["z"]["sy"] = chemin_son+"sondes_S_Y_CELL_VELOCITY_Z.son"
centered_velocity["z"]["sz"] = chemin_son+"sondes_S_Z_CELL_VELOCITY_Z.son"

indicatrice["sx"] = chemin_son+"sondes_S_X_INDICATRICE.son"
indicatrice["sy"] = chemin_son+"sondes_S_Y_INDICATRICE.son"
indicatrice["sz"] = chemin_son+"sondes_S_Z_INDICATRICE.son"


# Coordonnées + Plots
color=['r','g','b']

fig_sx, (ax_sx) = plt.subplots(1,1,figsize=(10,10))
fig_sy, (ax_sy) = plt.subplots(1,1,figsize=(10,10))
fig_sz, (ax_sz) = plt.subplots(1,1,figsize=(10,10))

fig_sx.suptitle("Sonde X - Plan XY")
fig_sy.suptitle("Sonde Y - Plan XY")
fig_sz.suptitle("Sonde Z - Plan XZ")
for i_vel, dir_vel in enumerate(["x","y","z"]):
    for i_sonde, dir_sonde in enumerate(["sx","sy","sz"]):
        coord_velocity[dir_vel][dir_sonde] = coordSonde(velocity[dir_vel][dir_sonde])
        coord_centered_velocity[dir_vel][dir_sonde] = coordSonde(centered_velocity[dir_vel][dir_sonde])
        coord_indicatrice[dir_sonde] = coordSonde(indicatrice[dir_sonde])
        
    ax_sx.scatter(         coord_velocity[dir_vel]["sx"]["x"],         coord_velocity[dir_vel]["sx"]["y"],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. coord.')
    ax_sx.scatter(coord_centered_velocity[dir_vel]["sx"]["x"],coord_centered_velocity[dir_vel]["sx"]["y"],color=color[i_vel],marker="x",label='cell_vel'+r'$_%s$'%(dir_vel)+'. coord.')
    ax_sx.scatter(               coord_indicatrice["sx"]["x"],               coord_indicatrice["sx"]["y"],color=color[i_vel],marker=".",label='ind. coord.')
                                                                                                                                       
    ax_sy.scatter(         coord_velocity[dir_vel]["sy"]["x"],         coord_velocity[dir_vel]["sy"]["y"],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. coord.')
    ax_sy.scatter(coord_centered_velocity[dir_vel]["sy"]["x"],coord_centered_velocity[dir_vel]["sy"]["y"],color=color[i_vel],marker="x",label='cell_vel'+r'_$%s$'%(dir_vel)+'. coord.')
    ax_sy.scatter(               coord_indicatrice["sy"]["x"],               coord_indicatrice["sy"]["y"],color=color[i_vel],marker=".",label='ind. coord.')
                                                                                                                                      
    ax_sz.scatter(         coord_velocity[dir_vel]["sz"]["x"],         coord_velocity[dir_vel]["sz"]["z"],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. coord.')
    ax_sz.scatter(coord_centered_velocity[dir_vel]["sz"]["x"],coord_centered_velocity[dir_vel]["sz"]["z"],color=color[i_vel],marker="x",label='cell_vel'+r'_$%s$'%(dir_vel)+'. coord.')
    ax_sz.scatter(               coord_indicatrice["sz"]["x"],               coord_indicatrice["sz"]["z"],color=color[i_vel],marker=".",label='ind. coord.')

ax_sx.set_xlabel('X')
ax_sy.set_xlabel('Y')
ax_sz.set_xlabel('X')
ax_sx.set_ylabel('Y')
ax_sy.set_ylabel('Y')
ax_sz.set_ylabel('Z')
ax_sx.grid(True,'both')
ax_sy.grid(True,'both')
ax_sz.grid(True,'both')
ax_sx.legend()
ax_sy.legend()
ax_sz.legend()

fig_sx.savefig("coord_sx_XY.png")
fig_sy.savefig("coord_sy_XY.png")
fig_sz.savefig("coord_sz_XZ.png")


# Valeurs relevées par les sondes

# Ajustements à enlever quand les sondes auront autant de points que de coordonnées
#  ou bien quand on aura compris pourquoi coord et valeur n'ont pas le mm nombre de valeurs...
# une fois ce contre temps regle, il faudra changer dans les scatter : cvxx par coor_ve...
cvxx  = {"x":0,"y":0,"z":0}
cvyy  = {"x":0,"y":0,"z":0}
cvzz  = {"x":0,"y":0,"z":0}
ccvxx  = {"x":0,"y":0,"z":0}
ccvyy  = {"x":0,"y":0,"z":0}
ccvzz  = {"x":0,"y":0,"z":0}
color=['r','g','b']

fig_sx, (ax_sx) = plt.subplots(1,1,figsize=(10,10))
fig_sy, (ax_sy) = plt.subplots(1,1,figsize=(10,10))
fig_sz, (ax_sz) = plt.subplots(1,1,figsize=(10,10))

fig_sx.suptitle("Valeurs Sonde X - Plan XY")
fig_sy.suptitle("Valeurs Sonde Y - Plan XY")
fig_sz.suptitle("Valeurs Sonde Z - Plan XZ")
for i_vel, dir_vel in enumerate(["x","y","z"]):
    for i_sonde, dir_sonde in enumerate(["sx","sy","sz"]):
        valeur_velocity[dir_vel][dir_sonde] =          loadPartSondeSeg(velocity[dir_vel][dir_sonde]         , 1,2)
        print(valeur_velocity[dir_vel][dir_sonde].shape)                                                           
        print(valeur_velocity[dir_vel][dir_sonde][0].shape) 
        print(coord_velocity[dir_vel]["sx"]["x"].shape)                                                          
        valeur_centered_velocity[dir_vel][dir_sonde] = loadPartSondeSeg(centered_velocity[dir_vel][dir_sonde], 1,2)
        print(valeur_centered_velocity[dir_vel][dir_sonde].shape)      
        print(valeur_centered_velocity[dir_vel][dir_sonde][0].shape) 
        print(coord_centered_velocity[dir_vel]["sx"]["x"].shape)       
        valeur_indicatrice[dir_sonde] =                loadPartSondeSeg(indicatrice[dir_sonde]               , 1,2)    
    
    cvxx[dir_vel]  = np.append(np.array(coord_velocity[dir_vel]["sx"]["x"]), 1.2*(coord_velocity[dir_vel]["sx"]["x"][len(coord_velocity[dir_vel]["sx"]["x"])-1])    )
    cvyy[dir_vel]  = np.append(np.array(coord_velocity[dir_vel]["sy"]["y"]), 1.2*(coord_velocity[dir_vel]["sy"]["y"][len(coord_velocity[dir_vel]["sy"]["y"])-1])    )
    cvzz[dir_vel]  = np.append(np.array(coord_velocity[dir_vel]["sz"]["z"]), 1.2*(coord_velocity[dir_vel]["sz"]["z"][len(coord_velocity[dir_vel]["sz"]["z"])-1])    )
    # cvxx[dir_vel]  = np.append(np.array(coord_velocity[dir_vel]["sx"]["x"]), (coord_velocity[dir_vel]["sx"]["x"][len(coord_velocity[dir_vel]["sx"]["x"])-1])    )
    # cvxx[dir_vel]  = np.append(np.array(coord_velocity[dir_vel]["sy"]["y"]), (coord_velocity[dir_vel]["sy"]["y"][len(coord_velocity[dir_vel]["sy"]["y"])-1])    )
    # cvxx[dir_vel]  = np.append(np.array(coord_velocity[dir_vel]["sz"]["z"]), (coord_velocity[dir_vel]["sz"]["z"][len(coord_velocity[dir_vel]["sz"]["z"])-1])    )
    ccvxx[dir_vel]  = np.append(np.array(coord_centered_velocity[dir_vel]["sx"]["x"]), 1.2*(coord_centered_velocity[dir_vel]["sx"]["x"][len(coord_centered_velocity[dir_vel]["sx"]["x"])-1])    )
    ccvyy[dir_vel]  = np.append(np.array(coord_centered_velocity[dir_vel]["sy"]["y"]), 1.2*(coord_centered_velocity[dir_vel]["sy"]["y"][len(coord_centered_velocity[dir_vel]["sy"]["y"])-1])    )
    ccvzz[dir_vel]  = np.append(np.array(coord_centered_velocity[dir_vel]["sz"]["z"]), 1.2*(coord_centered_velocity[dir_vel]["sz"]["z"][len(coord_centered_velocity[dir_vel]["sz"]["z"])-1])    )
    # ccvxx[dir_vel]  = np.append(np.array(coord_centered_velocity[dir_vel]["sx"]["x"]), (coord_centered_velocity[dir_vel]["sx"]["x"][len(coord_centered_velocity[dir_vel]["sx"]["x"])-1])    )
    # ccvxx[dir_vel]  = np.append(np.array(coord_centered_velocity[dir_vel]["sy"]["y"]), (coord_centered_velocity[dir_vel]["sy"]["y"][len(coord_centered_velocity[dir_vel]["sy"]["y"])-1])    )
    # ccvxx[dir_vel]  = np.append(np.array(coord_centered_velocity[dir_vel]["sz"]["z"]), (coord_centered_velocity[dir_vel]["sz"]["z"][len(coord_centered_velocity[dir_vel]["sz"]["z"])-1])    )
    # ax_sx.scatter(         coord_velocity[dir_vel]["sx"]["x"],         valeur_velocity[dir_vel]["sx"][0],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. valeur.')
    ax_sx.scatter(         cvxx[dir_vel],         valeur_velocity[dir_vel]["sx"][0],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. valeur.')
    ax_sx.scatter(ccvxx[dir_vel],valeur_centered_velocity[dir_vel]["sx"][0],color=color[i_vel],marker="x",label='cell_vel'+r'$_%s$'%(dir_vel)+'. valeur.')
                                                                                                                                       
    ax_sy.scatter(         cvyy[dir_vel],         valeur_velocity[dir_vel]["sy"][0],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. valeur.')
    ax_sy.scatter(ccvyy[dir_vel],valeur_centered_velocity[dir_vel]["sy"][0],color=color[i_vel],marker="x",label='cell_vel'+r'_$%s$'%(dir_vel)+'. valeur.')
                                                                                                                                     
    ax_sz.scatter(         cvzz[dir_vel],         valeur_velocity[dir_vel]["sz"][0],color=color[i_vel],marker="+",label='vel'+r'$_%s$'%(dir_vel)+'. valeur.')
    ax_sz.scatter(ccvzz[dir_vel],valeur_centered_velocity[dir_vel]["sz"][0],color=color[i_vel],marker="x",label='cell_vel'+r'_$%s$'%(dir_vel)+'. valeur.')

ax_sx.set_xlabel('X')
ax_sy.set_xlabel('Y')
ax_sz.set_xlabel('Z')
ax_sx.set_ylabel('Vel')
ax_sy.set_ylabel('Vel')
ax_sz.set_ylabel('Vel')
ax_sx.grid(True,'both')
ax_sy.grid(True,'both')
ax_sz.grid(True,'both')
ax_sx.legend()
ax_sy.legend()
ax_sz.legend()

fig_sx.savefig("valeur_sx_XY.png")
fig_sy.savefig("valeur_sy_XY.png")
fig_sz.savefig("valeur_sz_XZ.png")
plt.show()
