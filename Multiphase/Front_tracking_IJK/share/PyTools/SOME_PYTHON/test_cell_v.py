# -*- coding: utf-8 -*-                                                                            
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import sys
import os
import DNSTools3 as dtool

import math as mt
import time as time
import glob
import py_eav_II as pe
#import py_extract_all_volume

jdd_name = "cube"
repr_file = "NEXT/"+jdd_name

# Lecture vitesse relative
f=open("accplot.logF","r")
lines = f.readlines()
f.close()
for i, line in enumerate(lines):
    if "ur = " in line:
        ur = float(line.replace("ur = ",""))
# fin lecture vitesse relative

### DIRECTION X ###
print("### DIRECTION X ###")
Vx = pe.Champ()
liste_noms = glob.glob("VELOCITY_X_ELEM_DOM*.curve")
Vx.prepare_myself(liste_noms,jdd_name,repr_file)
# Chargement des champs
t = time.time()
Vx.load_fields()
t=time.time()-t
# fin chargement des champs
print("temps pour load_fields : ",t)
# Chargement des boites
t = time.time()
Vx.get_bubble_centers()
Vx.get_bubble_boxes()
Vx.get_field_around_bubble()
t=time.time()-t
# fin chargement des boites
print("temps pour chargement des  boites : ",t)

### DIRECTION Y ###
print("### DIRECTION Y ###")
Vy = pe.Champ()
liste_noms = glob.glob("VELOCITY_Y_ELEM_DOM*.curve")
Vy.prepare_myself(liste_noms,jdd_name,repr_file)
# Chargement des champs
t = time.time()
Vy.load_fields()
t=time.time()-t
# fin chargement des champs
print("temps pour load_fields : ",t)
# Chargement des boites
t = time.time()
Vy.get_bubble_centers()
Vy.get_bubble_boxes()
Vy.get_field_around_bubble()
t=time.time()-t
# fin chargement des boites
print("temps pour chargement des  boites : ",t)

### DIRECTION Z ###
print("### DIRECTION Z ###")
Vz = pe.Champ()
liste_noms = glob.glob("VELOCITY_Z_ELEM_DOM*.curve")
Vz.prepare_myself(liste_noms,jdd_name,repr_file)
# Chargement des champs
t = time.time()
Vz.load_fields()
t=time.time()-t
# fin chargement des champs
print("temps pour load_fields : ",t)
# Chargement des boites
t = time.time()
Vz.get_bubble_centers()
Vz.get_bubble_boxes()
Vz.get_field_around_bubble()
t=time.time()-t
# fin chargement des boites
print("temps pour chargement des  boites : ",t)


### PRESSION ###
print("### PRESSION ###")
Pr =  pe.Champ()
liste_noms = glob.glob("PRESSURE_ELEM_DOM*.curve")
Pr.prepare_myself(liste_noms,jdd_name,repr_file)
# Chargement des champs
t = time.time()
Pr.load_fields()
X_P,Y_P,Z_P = pe.Champ(),pe.Champ(),pe.Champ()
X_P.prepare_myself(liste_noms,jdd_name,repr_file)
Y_P.prepare_myself(liste_noms,jdd_name,repr_file)
Z_P.prepare_myself(liste_noms,jdd_name,repr_file)
[X_P.liste_field_kji, Y_P.liste_field_kji, Z_P.liste_field_kji] = Pr.elem_field_to_face_field()
t=time.time()-t
# fin chargement des champs
print("temps pour load_fields : ",t)
# Chargement des boites
t = time.time()
X_P.get_bubble_centers()
Y_P.get_bubble_centers()
Z_P.get_bubble_centers()
X_P.get_bubble_boxes()
Y_P.get_bubble_boxes()
Z_P.get_bubble_boxes()
X_P.get_field_around_bubble()
Y_P.get_field_around_bubble()
Z_P.get_field_around_bubble()
t=time.time()-t
# fin chargement des boites
print("temps pour chargement des  boites : ",t)


### INDICATRICE ###
print("### INDICATRICE ###")
Idc =  pe.Champ()
liste_noms = glob.glob("PRESSURE_ELEM_DOM*.curve")
Idc.prepare_myself(liste_noms,jdd_name,repr_file)
# Chargement des champs
t = time.time()
Idc.load_fields()
t=time.time()-t
# fin chargement des champs
print("temps pour load_fields : ",t)
### A REDIGER #####
# Correction fancy indicatrice on met a 100 uniquemen sur les cellules ne vlant ni 1 ni 0
# Idc.liste_data = [d+1 for d in liste_data]
# Idc.liste_field_kji = [f+1 for f in liste_field_kji]
###################
# Chargement des boites
t = time.time()
Idc.get_bubble_centers()
Idc.get_bubble_boxes()
Idc.get_field_around_bubble()
t=time.time()-t
# fin chargement des boites
print("temps pour chargement des  boites : ",t)

########################################################################
########################################################################
""" 
L'extraction est terminee, on passe au travail serieux.
Pour alleger le travail, on va certaineemnt redefinir des champs lights
pour lesquels on considere juste des plans de 5 mailles d'epaisseur,
autour du centre de la bulle.
"""
########################################################################
# Au cas ou, on a des fields allégés
# liste_bubble_field_tot_plan_kj = V.liste_bubble_field_kji[...,V.trailing_in_index-5:V.trailing_in_index+5]
# liste_bubble_field_tot_plan_ki = V.liste_bubble_field_kji[...,V.left_hand_in_index-5:V.left_hand_in_index+5,:]
# liste_bubble_field_tot_plan_ij = V.liste_bubble_field_kji[...,V.left_hand_in_index-5:V.left_hand_in_index+5]

# plt.figure(1)
# plt.imshow(liste_bubble_field_tot_plan_kj[4,0,:,:,5])
# plt.savefig("liste_bubble_field_tot_plan_kj.png")
# plt.figure(1)
# plt.imshow(liste_bubble_field_tot_plan_ki[4,0,:,5,:])
# plt.savefig("liste_bubble_field_tot_plan_ki.png")
# plt.figure(1)
# plt.imshow(liste_bubble_field_tot_plan_ij[4,0,5,:,:])
# plt.savefig("liste_bubble_field_tot_plan_ij.png")
########################################################################
V = [Vx,Vy,Vz]
P = [X_P, Y_P, Z_P]
nt,nb = Vx.liste_bubble_field_kji.shape[0],Vx.liste_bubble_field_kji.shape[1]
ncx,ncy,ncz = Vx.liste_bubble_field_kji.shape[2],Vx.liste_bubble_field_kji.shape[3],Vx.liste_bubble_field_kji.shape[4]

liste_bubble_tot_field_kji = np.zeros((3,nt,nb,ncx,ncy,ncz))
liste_bubble_bm_field_kji = np.zeros((3,nt,ncx,ncy,ncz))
liste_bubble_bm_tm_field_kji = np.zeros((3,ncx,ncy,ncz))
liste_bubble_pwf_field_kji = np.zeros((3,ncx,ncy,ncz))
liste_bubble_bia_field_kji = np.zeros((3,nt,nb,ncx,ncy,ncz))
liste_bubble_aib_field_kji = np.zeros((3,nt,nb,ncx,ncy,ncz))
liste_bubble_bia_p_times_v_field_kji = np.zeros((3,nt,nb,ncx,ncy,ncz))
########################################################################
t=time.time()
########################################################################
###### TRAVAIL - DOMAINE PHYSIQUE
for d,v in enumerate(V):
    ### CHAMPS DE VITESSE
    # total
    liste_bubble_tot_field_kji[d,...] = v.liste_bubble_field_kji.copy()+ ur
    # pwf
    # moyenne sur chaque bulle
    liste_bubble_bm_field_kji[d,...] = np.mean(v.liste_bubble_field_kji,axis=1)
    # ... et sur chaque pas de temps
    liste_bubble_bm_tm_field_kji[d,...] = np.mean(v.liste_bubble_field_kji,axis=0)
    liste_bubble_pwf_field_kji[d,...] = liste_bubble_bm_tm_field_kji[d,...]
    # bia
    liste_bubble_bia_field_kji[d,...] = liste_bubble_tot_field_kji[d,...] - liste_bubble_bm_tm_field_kji[d,...]
    liste_bubble_aib_field_kji[2-d,...] = liste_bubble_tot_field_kji[2-d,...] - liste_bubble_bm_tm_field_kji[2-d,...]
    
    ### CHAMP DE PRESSION*VITESSE
    liste_bubble_bia_p_times_v_field_kji[d,...] = P[d].liste_bubble_field_kji *liste_bubble_bia_field_kji[d,...]
    
### CHAMPS DE TENSEUR
# Rij BIA
liste_bubble_Rxx_field_kji = liste_bubble_bia_field_kji[0,...]**2
liste_bubble_Ryy_field_kji = liste_bubble_bia_field_kji[1,...]**2
liste_bubble_Rzz_field_kji = liste_bubble_bia_field_kji[2,...]**2
liste_bubble_Rxy_field_kji = liste_bubble_bia_field_kji[0,...]*liste_bubble_bia_field_kji[1,...]
liste_bubble_Rxz_field_kji = liste_bubble_bia_field_kji[0,...]*liste_bubble_bia_field_kji[2,...]
liste_bubble_Ryz_field_kji = liste_bubble_bia_field_kji[1,...]*liste_bubble_bia_field_kji[2,...]
# gradient  [direction_derive_kji,direction_vitesse,temps,bulle,k,j,i] :  d_kji v_ijk
if 'dual' in Vx.liste_nom_fichier:
    liste_bubble_grad_bia_field_kji = np.array(np.gradient(liste_bubble_bia_field_kji,Vx.L/(2*Vx.N),axis=(-3,-2,-1)))
else :
    liste_bubble_grad_bia_field_kji = np.array(np.gradient(liste_bubble_bia_field_kji,Vx.L/(Vx.N),axis=(-3,-2,-1)))
    
# gradient transpose
liste_bubble_grad_transpose_bia_field_kji = np.array(np.gradient(liste_bubble_aib_field_kji,Vx.L/(2*Vx.N),axis=(-1,-2,-3)))
# taux de deformations
liste_bubble_sij_bia_field_kji = 0.5*(liste_bubble_grad_bia_field_kji + liste_bubble_grad_transpose_bia_field_kji)
# taux de rotation
liste_bubble_oij_bia_field_kji = 0.5*(liste_bubble_grad_bia_field_kji - liste_bubble_grad_transpose_bia_field_kji)

### BILAN D'EC
# energie
liste_bubble_K_field_kji = 0.5 * (liste_bubble_Rxx_field_kji + liste_bubble_Ryy_field_kji + liste_bubble_Ryy_field_kji)
terme_K = np.mean(np.mean(liste_bubble_K_field_kji,axis=0),axis=0)
# Dissipation : 2 nu sij sij
nu = 0.00035
liste_bubble_dissipation = 2*nu*np.einsum("ij...,ij...->...",liste_bubble_sij_bia_field_kji,liste_bubble_sij_bia_field_kji)
terme_dissipation = np.mean(np.mean(liste_bubble_dissipation,axis=0),axis=0)
# Production = Rij * dj ui
liste_bubble_production = (liste_bubble_Rxx_field_kji*liste_bubble_grad_bia_field_kji[2,2,...] +
                           liste_bubble_Rxy_field_kji*liste_bubble_grad_bia_field_kji[1,0,...] + 
                           liste_bubble_Rxz_field_kji*liste_bubble_grad_bia_field_kji[2,0,...] + 
                           liste_bubble_Rxy_field_kji*liste_bubble_grad_bia_field_kji[0,1,...] + 
                           liste_bubble_Ryy_field_kji*liste_bubble_grad_bia_field_kji[1,1,...] + 
                           liste_bubble_Ryz_field_kji*liste_bubble_grad_bia_field_kji[1,2,...] + 
                           liste_bubble_Rxz_field_kji*liste_bubble_grad_bia_field_kji[2,0,...] + 
                           liste_bubble_Ryz_field_kji*liste_bubble_grad_bia_field_kji[2,1,...] + 
                           liste_bubble_Rzz_field_kji*liste_bubble_grad_bia_field_kji[0,0,...] )
terme_production = np.mean(np.mean(liste_bubble_production,axis=0),axis=0)
# Transfert total d'énergie (pression-diffusion, convection, transfert non-lineaire )
liste_bubble_transfer = liste_bubble_production - liste_bubble_dissipation
terme_transfer = np.mean(np.mean(liste_bubble_transfer,axis=0),axis=0)
# Pression-diffusion : dj (p/rho uj)
liste_bubble_pression_diffusion = (np.gradient(liste_bubble_bia_p_times_v_field_kji[0,...],Vx.L/(Vx.N),axis=(-3,-2,-1))[0] +
                                   np.gradient(liste_bubble_bia_p_times_v_field_kji[1,...],Vx.L/(Vx.N),axis=(-3,-2,-1))[1] +
                                   np.gradient(liste_bubble_bia_p_times_v_field_kji[2,...],Vx.L/(Vx.N),axis=(-3,-2,-1))[2] )
terme_pression_diffusion = np.mean(np.mean(liste_bubble_pression_diffusion,axis=0),axis=0)
# Transfer non-linéaire : 1/2 dj (ui ui uj)
liste_bubble_transfert_non_lineaire = (np.gradient(liste_bubble_bia_field_kji[0,...]**3,Vx.L/(Vx.N),axis=(-1)) +
                                       np.gradient(liste_bubble_bia_field_kji[0,...]**2 * liste_bubble_bia_field_kji[1,...],Vx.L/(Vx.N),axis=(-2)) +
                                       np.gradient(liste_bubble_bia_field_kji[0,...]**2 * liste_bubble_bia_field_kji[2,...],Vx.L/(Vx.N),axis=(-3)) +
                                       np.gradient(liste_bubble_bia_field_kji[1,...]**3,Vx.L/(Vx.N),axis=(-2)) +
                                       np.gradient(liste_bubble_bia_field_kji[1,...]**2 * liste_bubble_bia_field_kji[0,...],Vx.L/(Vx.N),axis=(-1)) +
                                       np.gradient(liste_bubble_bia_field_kji[1,...]**2 * liste_bubble_bia_field_kji[2,...],Vx.L/(Vx.N),axis=(-3)) +
                                       np.gradient(liste_bubble_bia_field_kji[2,...]**3,Vx.L/(Vx.N),axis=(-3)) +
                                       np.gradient(liste_bubble_bia_field_kji[2,...]**2 * liste_bubble_bia_field_kji[0,...],Vx.L/(Vx.N),axis=(-1)) +
                                       np.gradient(liste_bubble_bia_field_kji[2,...]**2 * liste_bubble_bia_field_kji[1,...],Vx.L/(Vx.N),axis=(-2)) )
terme_transfert_non_lineaire = np.mean(np.mean(liste_bubble_transfert_non_lineaire,axis=0),axis=0)
                                      
# Convection : dj(<uj> K)
liste_bubble_transfert_convectif = (np.gradient(liste_bubble_K_field_kji*liste_bubble_bia_field_kji[0,...],Vx.L/(Vx.N),axis=(-3,-2,-1))[0] +
                                    np.gradient(liste_bubble_K_field_kji*liste_bubble_bia_field_kji[1,...],Vx.L/(Vx.N),axis=(-3,-2,-1))[1] +
                                    np.gradient(liste_bubble_K_field_kji*liste_bubble_bia_field_kji[2,...],Vx.L/(Vx.N),axis=(-3,-2,-1))[2] )
terme_transfert_convectif = np.mean(np.mean(liste_bubble_transfert_convectif,axis=0),axis=0)


########################################################################
########################################################################
t=time.time()-t

########################################################################
###       GRANDEURS NON MOYENNEES, INSTANTANNES                      ###
########################################################################
indice_kj = 17
indice_i = 22
plt.figure(1)
""" [composante, iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_bia_field_kji[0,4,0,indice_kj,:,:])
plt.savefig("BIA_xy.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_K_field_kji[4,0,indice_kj,:,:])
plt.savefig("K_xy.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_dissipation[4,0,indice_kj,:,:])
plt.savefig("Dissip_xy.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_production[4,0,indice_kj,:,:])
plt.savefig("Production.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_transfer[4,0,indice_kj,:,:])
plt.savefig("Transfert_total_xy.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_pression_diffusion[4,0,indice_kj,:,:])
plt.savefig("Pression_diffusion_xy.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(liste_bubble_transfert_convectif[4,0,indice_kj,:,:])
plt.savefig("Convection_xy.png")

plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(Vx.liste_bubble_field_kji[4,0,indice_kj,:,:])
plt.savefig("Utot_xy.png")


########################################################################
###       GRANDEURS MOYENNEES                                        ###
########################################################################


plt.figure(1)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(terme_K[indice_kj,:,:], cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.savefig("mK_xy.png")

plt.figure(2)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(terme_dissipation[indice_kj,:,:], cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.savefig("mDissip_xy.png")

plt.figure(3)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(terme_production[indice_kj,:,:], cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.savefig("mProduction.png")

plt.figure(4)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(terme_transfer[indice_kj,:,:], cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.savefig("mTransfert_total_xy.png")

plt.figure(5)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(terme_pression_diffusion[indice_kj,:,:], cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.savefig("mPression_diffusion_xy.png")

plt.figure(6)
""" [iteration, bulle, cellule_k, cellule_j, cellule_i] """
plt.imshow(terme_transfert_convectif[indice_kj,:,:], cmap=plt.cm.BuPu_r)
plt.colorbar()
plt.savefig("mConvection_xy.png")






print("temps pour TRAVAIL : ",t)


"""
Il reste encore a reprendre ce que je faisais pour tous les termes du bilan
de tke dans l'espace physique
Il reste encore a finir d'ecrire-adapter ce que j'ai fais pour le bilan de
tke dans le domaine spectral
"""

