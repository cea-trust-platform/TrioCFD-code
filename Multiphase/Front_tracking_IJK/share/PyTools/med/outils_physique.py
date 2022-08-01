# -*- coding: utf-8 -*-
from optparse import OptionParser
import MEDCouplingRemapper # Pour l'interpolateur
import numpy as np
import time # Pour la création d'un nouveau maillage irrégulier
import MEDLoader as ml
import methodes_2020 # comme ca on est certain que si on prend juste une fonction
                     # qui en appelait d'autres de ce mm fichier, elle ocntinuera de ofnctionner.
import outils_physique
import matplotlib.pyplot as plt


def grad_scal(f,BC='perioX, perioY, perioZ'):
    # Permet de calculer le gradient d'un champ f à UNE composante
    # Dans le cas de la calibration de KT on calcule uniquement des gradients de température
    # f ne doit pas forcément être associé à un maillage cartésien.
    # f doit etre du meme type que ce que renvoie LoadFieldConservative
    # f peut etre une sortie de getFieldOnMeshAtLevel()
    # BC specifie le type de condition limite que l'on impose

    ####################################################################
    ################## TRAVAIL PREPARATOIRE ############################
    f2,mC = methodes_2020.en_tete(f)
    xn,yn,zn,champ = methodes_2020.MEDFieldToNumpy(f2,mC)
    Nx = len(xn) - 1;    Ny = len(yn) - 1;    Nz = len(zn) - 1
    grad = np.zeros((Nx, Ny, Nz, 3)) # Gradient = 3 compo...
    ####################################################################


    
    ####################################################################
    ##################                     #############################
    ############        COEUR DU GRADIENT         ######################
    ##################                     #############################
    ####################################################################
    
    ####################################################################
    ############        GRAD COMPO X :            ######################
    ####################################################################
    
    # PATH : S'il n'y a pas assez de mailles en x... 
    if grad.shape[0]< 3: 
        print("Pas suffisament de maille selon x pour calculer le gradient")
        pass
    elif grad.shape[0]>= 3:
        # Calcul des pas du maillages selon x :
        v = np.array(xn.getValues()) # Conversion en NumPy array rappel : xn = mC.getCoordsAt(0) (les noeuds)
        dx = v[1:] - v[:-1] # dx_i = x_i+1 - x_i
        # L'evaluation du dx entre 2 cellules :
        h = (dx[:-1] + dx[1:]) / 2.# contient (Nsommets - 2) valeurs = Ncellules - 1


        if not('perioX' in BC):
            ################################################################
            # Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[:,np.newaxis,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[:,np.newaxis,np.newaxis]        
            u_pl = champ[2:,  :,:,0]
            u_m = champ[:-2,  :,:,0]
            u_c = champ[1:-1, :,:,0]
            grad[1:-1,:,:, 0] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord gauche :
            ################################################################
            h1 = h[0] 
            h2 = h[1]
            u_g   = champ[0,:,:,0]
            u_pl  = champ[1,:,:,0]
            u_pl2 = champ[2,:,:,0]
            grad[0,:,:,0] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord droit : 
            ################################################################
            h1 = h[-2] 
            h2 = h[-1]
            u_d  = champ[-1,:,:,0]
            u_m  = champ[-2,:,:,0]
            u_m2 = champ[-3,:,:,0]
            grad[-1,:,:,0] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
            pass
        else :
            ################################################################
            # PERIODIQUE : Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h = np.concatenate((h[-1],h),axis=None)
            h = np.concatenate((h,h[1]),axis=None)
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[:,np.newaxis,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[:,np.newaxis,np.newaxis]

            bout_u_pl = np.reshape([champ[0,:,:,0]],(1,Ny,Nz))
            bout_u_m =  np.reshape([champ[-1,:,:,0]],(1,Ny,Nz))
            
            u_pl = np.concatenate((champ[1:,:,:,0], bout_u_pl),axis=0)
            u_m = np.concatenate((bout_u_m ,champ[:-1,  :,:,0]),axis=0)
            u_c = champ[:, :,:,0]
            grad[:,:,:, 0] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)

            
    ####################################################################
    ############        GRAD COMPO Y :            ######################
    ####################################################################

    # PATH : S'il n'y a pas assez de mailles en y...
    if grad.shape[1]< 3: 
        print("Pas suffisament de maille selon y pour calculer le gradient")
        pass
    elif grad.shape[1]>= 3:
        # Calcul des pas du maillages selon y :
        v = np.array(yn.getValues()) # Conversion en NumPy array
        dy = v[1:] - v[:-1] # dy_i = y_i+1 - y_i
        # L'evaluation du dx entre 2 cellules :
        h = (dy[:-1] + dy[1:]) / 2.# contient (Nsommets - 2) valeurs = Ncellules - 1


        if not ('perioY' in BC):
            ################################################################
            # Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################        
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[:,np.newaxis,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[:,np.newaxis,np.newaxis]
            u_pl = champ[:,2:,:,0]
            u_m = champ[:,:-2,:,0]
            u_c = champ[:,1:-1,:,0]
            grad[:,1:-1,:, 1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord gauche :
            ################################################################
            h1 = h[0] 
            h2 = h[1]
            u_g   = champ[:,0,:,0]
            u_pl  = champ[:,1,:,0]
            u_pl2 = champ[:,2,:,0]
            grad[:,0,:,1] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord droit : 
            ################################################################
            h1 = h[-2] 
            h2 = h[-1]
            u_d  = champ[:,-1,:,0]
            u_m  = champ[:,-2,:,0]
            u_m2 = champ[:,-3,:,0]
            grad[:,-1,:,1] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
            pass
        else :
            ################################################################
            # PERIODIQUE : Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h = np.concatenate((h[-1],h),axis=None)
            h = np.concatenate((h,h[1]),axis=None)
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[np.newaxis,:,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[np.newaxis,:,np.newaxis]

            bout_u_pl = np.reshape([champ[:,0,:,0]],(Nx,1,Nz))
            bout_u_m = np.reshape([champ[:,-1,:,0]],(Nx,1,Nz))
        
            u_pl = np.concatenate((champ[:,1:,:,0], bout_u_pl),axis=1)
            u_m = np.concatenate((bout_u_m,champ[:,:-1,:,0]), axis=1)
            u_c = champ[:, :,:,0]
            grad[:,:,:, 1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            
            
    ####################################################################
    ############        GRAD COMPO Z :            ######################
    ####################################################################

    # PATH : S'il n'y a pas assez de mailles en z...
    if grad.shape[2]< 3: 
        print("Pas suffisament de maille selon z pour calculer le gradient")
        pass
    elif grad.shape[2]>= 3:
        # Calcul des pas du maillages selon z :
        v = np.array(zn.getValues()) # Conversion en NumPy array
        dz = v[1:] - v[:-1] # dz_i = z_i+1 - z_i
        # L'evaluation du dx entre 2 cellules :
        h = (dz[:-1] + dz[1:]) / 2. # contient (Nsommets - 2) valeurs = Ncellules - 1


        if not ('perioZ' in BC):
            ################################################################
            # Formule centrée (ordre 2) pour pas variable dans le domaine :
            ################################################################
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[np.newaxis,np.newaxis, :] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[np.newaxis,np.newaxis, :]
            u_pl = champ[:,:,2:,  0]
            u_m = champ[:,:,:-2,  0]
            u_c = champ[ :,:,1:-1,0]
            grad[... ,1:-1, 2] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord gauche :
            ################################################################
            h1 = h[0] 
            h2 = h[1]
            u_g   = champ[:,:,0,0]
            u_pl  = champ[:,:,1,0]
            u_pl2 = champ[:,:,2,0]
            grad[:,:,0,2] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord droit : 
            ################################################################
            h1 = h[-2] 
            h2 = h[-1]
            u_d  = champ[:,:,-1,0]
            u_m  = champ[:,:,-2,0]
            u_m2 = champ[:,:,-3,0]
            grad[:,:,-1,2] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
            pass
        else :
            ################################################################
            # PERIODIQUE : Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h = np.concatenate((h[-1],h),axis=None)
            h = np.concatenate((h,h[1]),axis=None)
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[np.newaxis,np.newaxis,:] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[np.newaxis,np.newaxis,:]

            bout_u_pl = np.reshape([champ[:,:,0,0]],(Nx,Ny,1))
            bout_u_m =  np.reshape([champ[:,:,-1,0]],(Nx,Ny,1))
            
            u_pl = np.concatenate((champ[:,:,1:,0],bout_u_pl), axis = 2)
            u_m =  np.concatenate((bout_u_m,champ[:,:,:-1,0]),axis = 2)
            u_c = champ[:,:, :,0]
          
            grad[:,:,:, 2] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            pass

    ####################################################################
    # Passage du NumPy à MED
    fgrad = methodes_2020.NumpyToMEDField(grad, mC, NomChamp="Gradient")
    ####################################################################
    
    return grad, fgrad


def grad_vect(f,composante,BC='perioX, perioY, perioZ'):

    ####################################################################
    ################## TRAVAIL PREPARATOIRE ############################
    f2,mC = methodes_2020.en_tete(f)
    xn,yn,zn,champ = methodes_2020.MEDFieldToNumpy(f2,mC)
    Nx = len(xn) - 1;    Ny = len(yn) - 1;    Nz = len(zn) - 1
    grad = np.zeros((Nx, Ny, Nz, 3)) # Gradient = 3 compo...
    ####################################################################


    
    ####################################################################
    ##################                     #############################
    ############        COEUR DU GRADIENT         ######################
    ##################                     #############################
    ####################################################################
    
    ####################################################################
    ############        GRAD COMPO X :            ######################
    ####################################################################
    
    # PATH : S'il n'y a pas assez de mailles en x...
    if grad.shape[0]< 3: 
        print("Pas suffisament de maille selon x pour calculer le gradient")
        pass
    elif grad.shape[0]>= 3:
        # Calcul des pas du maillages selon x :   
        v = np.array(xn.getValues()) # Conversion en NumPy array
        dx = v[1:] - v[:-1] # dx_i = x_i+1 - x_i
        # L'evaluation du dx entre 2 cellules :
        h = (dx[:-1] + dx[1:]) / 2.# contient (Nsommets - 2) valeurs = Ncellules - 1

        
        if not('perioX' in BC):
            ################################################################
            # Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[:,np.newaxis,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[:,np.newaxis,np.newaxis]        
            u_pl = champ[2:,  :,:,composante]
            u_m = champ[:-2,  :,:,composante]
            u_c = champ[1:-1, :,:,composante]
            grad[1:-1,:,:, 0] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord gauche :
            ################################################################
            h1 = h[0] 
            h2 = h[1]
            u_g   = champ[0,:,:,composante]
            u_pl  = champ[1,:,:,composante]
            u_pl2 = champ[2,:,:,composante]
            grad[0,:,:,0] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord droit : 
            ################################################################
            h1 = h[-2] 
            h2 = h[-1]
            u_d  = champ[-1,:,:,composante0]
            u_m  = champ[-2,:,:,composante0]
            u_m2 = champ[-3,:,:,composante0]
            grad[-1,:,:,0] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
            pass
        else :
            ################################################################
            # PERIODIQUE : Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h = np.concatenate((h[-1],h),axis=None)
            h = np.concatenate((h,h[1]),axis=None)
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[:,np.newaxis,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[:,np.newaxis,np.newaxis]

            bout_u_pl = np.reshape([champ[0,:,:,composante]],(1,Ny,Nz))
            bout_u_m =  np.reshape([champ[-1,:,:,composante]],(1,Ny,Nz))
            
            u_pl = np.concatenate((champ[1:,:,:,composante], bout_u_pl),axis=0)
            u_m = np.concatenate((bout_u_m ,champ[:-1,  :,:,composante]),axis=0)
            u_c = champ[:, :,:,composante]
            grad[:,:,:, 0] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)


    ####################################################################
    ############        GRAD COMPO Y :            ######################
    ####################################################################

    # PATH : S'il n'y a pas assez de mailles en y...
    if grad.shape[1]< 3: 
        print("Pas suffisament de maille selon y pour calculer le gradient")
        pass
    elif grad.shape[1]>= 3:
        # Calcul des pas du maillages selon y :
        v = np.array(yn.getValues()) # Conversion en NumPy array
        dy = v[1:] - v[:-1] # dy_i = y_i+1 - y_i
        # L'evaluation du dx entre 2 cellules :
        h = (dy[:-1] + dy[1:]) / 2.# contient (Nsommets - 2) valeurs = Ncellules - 1


        if not ('perioY' in BC):
            ################################################################
            # Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################        
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[:,np.newaxis,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[:,np.newaxis,np.newaxis]
            u_pl = champ[:,2:,:,composante]
            u_m = champ[:,:-2,:,composante]
            u_c = champ[:,1:-1,:,composante]
            grad[:,1:-1,:, 1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord gauche :
            ################################################################
            h1 = h[0] 
            h2 = h[1]
            u_g   = champ[:,0,:,composante]
            u_pl  = champ[:,1,:,composante]
            u_pl2 = champ[:,2,:,composante]
            grad[:,0,:,1] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord droit : 
            ################################################################
            h1 = h[-2] 
            h2 = h[-1]
            u_d  = champ[:,-1,:,composante]
            u_m  = champ[:,-2,:,composante]
            u_m2 = champ[:,-3,:,composante]
            grad[:,-1,:,1] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
            pass
        else :
            ################################################################
            # PERIODIQUE : Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h = np.concatenate((h[-1],h),axis=None)
            h = np.concatenate((h,h[1]),axis=None)
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[np.newaxis,:,np.newaxis] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[np.newaxis,:,np.newaxis]

            bout_u_pl = np.reshape([champ[:,0,:,composante]],(Nx,1,Nz))
            bout_u_m = np.reshape([champ[:,-1,:,composante]],(Nx,1,Nz))

            u_pl = np.concatenate((champ[:,1:,:,composante], bout_u_pl),axis=1)
            u_m = np.concatenate((bout_u_m,champ[:,:-1,:,composante]), axis=1)
            u_c = champ[:, :,:,composante]

            grad[:,:,:, 1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)


            

    ####################################################################
    ############        GRAD COMPO Z :            ######################
    ####################################################################


    # PATH : S'il n'y a pas assez de mailles en z...
    if grad.shape[2]< 3: 
        print("Pas suffisament de maille selon z pour calculer le gradient")
        pass
    elif grad.shape[2]>= 3:
        # Calcul des pas du maillages selon z :
        v = np.array(zn.getValues()) # Conversion en NumPy array
        dz = v[1:] - v[:-1] # dz_i = z_i+1 - z_i
        # L'evaluation du dx entre 2 cellules :
        h = (dz[:-1] + dz[1:]) / 2. # contient (Nsommets - 2) valeurs = Ncellules - 1

        
        if not ('perioZ' in BC):
            ################################################################
            # Formule centrée (ordre 2) pour pas variable dans le domaine :
            ################################################################
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[np.newaxis,np.newaxis, :] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[np.newaxis,np.newaxis, :]
            u_pl = champ[:,:,2:,  composante]
            u_m = champ[:,:,:-2,  composante]
            u_c = champ[ :,:,1:-1,composante]
            grad[... ,1:-1, 2] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord gauche :
            ################################################################
            h1 = h[0] 
            h2 = h[1]
            u_g   = champ[:,:,0,composante]
            u_pl  = champ[:,:,1,composante]
            u_pl2 = champ[:,:,2,composante]
            grad[:,:,0,2] = (1/(h1*h2*(h1+h2))) * (-(2*h1*h2+h2**2)*u_g + ((h1+h2)**2)*u_pl - (h1**2)*u_pl2)
            ################################################################
            # Formule décentrée (ordre 2) pour pas variable sur le bord droit : 
            ################################################################
            h1 = h[-2] 
            h2 = h[-1]
            u_d  = champ[:,:,-1,composante]
            u_m  = champ[:,:,-2,composante]
            u_m2 = champ[:,:,-3,composante]
            grad[:,:,-1,2] = 1/(h1*h2*(h1+h2)) * ((h2**2)*u_m2 - ((h1+h2)**2)*u_m - ((h2**2)-(h1+h2)**2)*u_d)
            pass
        else :
            ################################################################
            # PERIODIQUE : Formule centrée (ordre 2) pour pas variable dans le domaine : 
            ################################################################
            h = np.concatenate((h[-1],h),axis=None)
            h = np.concatenate((h,h[1]),axis=None)
            h1 = h[:-1]
            h2 = h[1:]
            h1 = h1[np.newaxis,np.newaxis,:] # Passage d'un vecteur à une matrice 3D.
            h2 = h2[np.newaxis,np.newaxis,:]

            bout_u_pl = np.reshape([champ[:,:,0,composante]],(Nx,Ny,1))
            bout_u_m =  np.reshape([champ[:,:,-1,composante]],(Nx,Ny,1))
            
            u_pl = np.concatenate((champ[:,:,1:,composante],bout_u_pl), axis = 2)
            u_m =  np.concatenate((bout_u_m,champ[:,:,:-1,composante]),axis = 2)
            u_c = champ[:,:, :,composante]
           
            grad[:,:,:, 2] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
            pass
    return grad
    

def jacob(f,BC='perioX,perioY,perioZ'):
    ####################################################################
    # Permet de calculer la jacobienne d'un champ f à 3 composantes
    # f : champ MED à 3 composantes (pas plus, pas moins)
    # jacob  : jacobienne de f, format NumPy
    # fjacob : jacobienne de f, format MEDField
    ####################################################################
    
    # ###################################################################
    # ################# TRAVAIL PREPARATOIRE ############################
    f2,mC = methodes_2020.en_tete(f)
    xn,yn,zn = methodes_2020.MEDFieldToNumpy_SansChamp(mC)
    Nx = len(xn) - 1;    Ny = len(yn) - 1;    Nz = len(zn) - 1
    jacob = np.zeros((Nx, Ny, Nz, 3, 3)) # Jacobienne = 9 compo...
    # ###################################################################

    for i in range(3):
        composante = i
        gradient = outils_physique.grad_vect(f2, composante,BC)
        for l in range(Nx):
            for m in range(Ny):
                for n in range(Nz):
                    jacob[l,m,n,composante] = gradient[l,m,n]
                    pass
                pass
            pass

    # ###################################################################
    # Passage du NumPy à MED
    fjacob = methodes_2020.NumpyToMEDField(jacob,mC, NomChamp = 'grad_'+f.getName())
    # ###################################################################
    return jacob, fjacob

# Delta est la demi-taille du filtre : 
def filtre(champ, mccfd, delta, BC='perioX, perioY, perioZ'):
    # @champ : champ au format numpy (shape : (Nx,Ny,Nz,Ncompo)) ou au format MED
    # @mccfd : maillage cartésien sur lequel s'appuye le champ
    # @delta : demi-taille du fitre à appliquer

    if type(champ)!=type(np.array(351)):
        print("\n\n CE N'EST PAS UN ARRAY !!!")
        f2 = champ.deepCopy()
        f2,mC = methodes_2020.en_tete(f2)
        champ = f2.deepCopy()
        xn,yn,zn,champ = methodes_2020.MEDFieldToNumpy(champ,mC)
        nx = len(xn) - 1
        ny = len(yn) - 1
        nz = len(zn) - 1
        Nco = f2.getNumberOfComponents()


    else :
        mC = mccfd
        pass

    
    # Récupération des corrdonnées des noeuds du maillage
    noeudx = mccfd.getCoordsAt(0).toNumPyArray()
    noeudy = mccfd.getCoordsAt(1).toNumPyArray()
    noeudz = mccfd.getCoordsAt(2).toNumPyArray()
    #
    Nx = len(noeudx)-1
    Ny = len(noeudy)-1
    Nz = len(noeudz)-1
    # Récupération des coordonnées des centre des cellules
    baryx = (noeudx[1:] + noeudx[:-1])/2
    baryy = (noeudy[1:] + noeudy[:-1])/2
    baryz = (noeudz[1:] + noeudz[:-1])/2
    #
    # Calcul de la longeur des cotés des cellules
    dx = (noeudx[1:] - noeudx[:-1])
    dxx = dx.reshape(Nx,1,1)
    dy = (noeudy[1:] - noeudy[:-1])
    dyy = dy.reshape(1,Ny,1)
    dz = (noeudz[1:] - noeudz[:-1])
    dzz = dz.reshape(1,1,Nz)
    #
    # Calcul du volume des cellules
    vol = dxx*dyy*dzz
    #
    # Création des tableaux mins et maxs
    mini_x, maxi_x = methodes_2020.get_mini_maxi(baryx,delta)
    mini_y, maxi_y = methodes_2020.get_mini_maxi(baryy,delta)
    mini_z, maxi_z = methodes_2020.get_mini_maxi(baryz,delta)
    #print "x : "  mini_x - maxi_x
    #print "y : "  mini_y - maxi_y
    #print "z : "  mini_z - maxi_z
    #
    sh = np.shape(champ)
    print("shape du champ")
    print(sh)
    if (len(sh) == 5): 
        # Cas tenseur : 
        resu = np.zeros(sh)
        for compo in range(sh[4]):
            tmp = champ[:,:,:,:,compo]
            resu[:,:,:,:,compo] = methodes_2020.filtre_k(tmp, vol, mini_x, maxi_x, mini_y, maxi_y, mini_z, maxi_z)
            pass
        # Passage du NumPy à MED
        fresu = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
        fresu.setNature(ml.IntensiveMaximum)
        fresu.setName("filtre")
        fresu.setMesh(mccfd)
        fresu.fillFromAnalytic(sh[4]*sh[3], "0")
        trans = np.transpose(resu, (2,1,0,3,4))   # Nz, Ny, Nx, Ncomposante, NtermeParComposante
        flat = trans.flatten()
        arr = ml.DataArrayDouble(flat)
        arr.rearrange(sh[4]*sh[3])
        fresu.setArray(arr)
        #
    else : 
        #Cas scalaire ou vecteur :
        resu = methodes_2020.filtre_k(champ, vol, mini_x, maxi_x, mini_y, maxi_y, mini_z, maxi_z)

        # Passage du NumPy à MED
        fresu = ml.MEDCouplingFieldDouble(ml.ON_CELLS)
        fresu.setNature(ml.IntensiveMaximum)
        fresu.setName("filtre")
        fresu.setMesh(mccfd)
        fresu.fillFromAnalytic(sh[3], "0")
        trans = np.transpose(resu, (2,1,0,3))   # Nz, Ny, Nx, Ncomposante, NtermeParComposante
        flat = trans.flatten()
        arr = ml.DataArrayDouble(flat)
        arr.rearrange(sh[3])
        fresu.setArray(arr)
        ### AAAARG SOUCIS CAR PAS DE MESH : PTET mccfd  

    return resu, fresu


"""
    Bon, comme les trois foncitons précédents, sont vraiment pas évidentes à comprendre de fond
    en comble, qu'on a un résultat pas physique, mais qu'en fait ca vient de get_mini_maxi, qui
    pourtant fait son boulot comme on a compris qu'il devait le faire, et surtout que filtre
    contient une triple boucle imbriquée chacune sur le nombre de cellules du domaine , on va plutot refaire ce filtrage. Et si ça ne convient pas, he ben on essaiera plus fort de fair efoncitonner les trois fonction sprécédentes. Pour le moment nle nouveau filtrage va comprendre une triple boucle imbriquee sur le nombre de directions.
"""

def filtre_gab(delta,champ,tailles,BC='perioX,perioY,perioZ'):
    """
        CA C PAS FOU : CA PREND BEAUCOUP DE MEMOIRE EN FAIT
        QUAND MEME DEMANDER DE REGARDER UN PEU PCK JE SENS QUE CA PEUT
        ETRE PAS MAUVAIS EN VRAI
    """
    """
    delta : float --> demi-taille du filtre de moyennage
    champ : array de float --> le champ à traiter
    tailles : array [Lx, Ly, Lz] --> La longueur de la boite dan schacune des
    dimensions.

    On dit que le maillage est régulier, mais une maille peut etre non cubique

    POUR LE MOMENT ON SE PLACE EN 3D PERIODIQUE
    """
    sh = np.shape(champ)
    print "sh", sh
    # Passage d'une longueur delta à un nombre d'indices
    Nx, Ny, Nz = sh[0], sh[1], sh[2]
    Lx, Ly, Lz = tailles[0], tailles[1], tailles[2]
    dx, dy, dz = Lx/Nx, Ly/Ny, Lz/Nz
    
    shiftTotX = int(delta*Nx/Lx) + 1
    shiftTotY = int(delta*Ny/Ly) + 1
    shiftTotZ = int(delta*Nz/Lz) + 1
    # print("delta,Nx,Ny,Nz,Lx,Ly,Lz",delta,Nx,Ny,Nz,Lx,Ly,Lz)
    # print "shiftX,shiftY,shiftZ",shiftTotX,shiftTotY,shiftTotZ

    # Calcul du volume d'une maille. Maillage régulier donc toutes les mailles
    # ont le mm volume !
    vol = np.prod(tailles)/np.prod(sh[:-1])         # prod : comme sum mais avec la multiplication
    print("volume d'une cellule",vol)    
    
    champ_tout_shift = champ.copy() * vol
    Nchamp = 1
    shiftDirection = []    
    shiftX,shiftY,shiftZ = 0,0,0
    while shiftX < shiftTotX or shiftY < shiftTotY or shiftZ < shiftTotZ :
        if shiftX < shiftTotX :
            shiftX+=1
            if not ("X" in shiftDirection) :
                shiftDirection.append("X")
        if shiftY < shiftTotY :
            shiftY+=1
            if not ("Y" in shiftDirection) :
                shiftDirection.append("Y")
        if shiftZ < shiftTotZ :
            shiftZ+=1
            if not ("Z" in shiftDirection) :
                shiftDirection.append("Z")
        
        Shift = shiftX,shiftY,shiftZ
        # print('shiftDirection,Shift : ',shiftDirection,Shift)
        champ_un_shift,NUnShift = methodes_2020.filtre_un_shift(Shift,champ,shiftDirection,vol)
        # print('NUnShift',NUnShift)
        Nchamp += NUnShift
        champ_tout_shift += champ_un_shift

    # volume total du filtre
    Vol = vol * (Nchamp)
    # print "Volume de tout le filtre", Vol

    # Nombre de cellules dans le filtre
    # print("Nombre de cellules dans le filtre", Nchamp)

    # On moyenne le tout 
    # champ_tout_shift /= Vol

    return (champ_tout_shift/Vol)


def Applatit(MEDchamp,temps,direction):
    ####################################################################
    # MEDchamp    : champ MED a moyenner
    #               (Nx,Ny,Nz,Ncompo)
    # direction   : direction selon laquelle on moyenne
    #               pour direction = "X" : (1 ,Ny,Nz,Ncompo)
    # temps       : [temps, iteration, order]
    # --> pour avoir les iterations : MEDchamp.getIterations()
    ####################################################################
    
    # Preparation ######################################################
    time = MEDchamp.getTime()
    MEDchamp,mC    = methodes_2020.en_tete(MEDchamp)
    xn,yn,zn,champ = methodes_2020.MEDFieldToNumpy(MEDchamp,mC)

    # Coeur de fonction ################################################
    sh = list(champ.shape)
    ax = 25 #to raise an error if the direction is incorrectly given.
    ax = (0*(direction=="X") + 1*(direction=="Y") + 2*(direction=="Z"))
    sh[ax] = 1
    champPlat = champ.mean(axis=ax).reshape(sh)

    # MED vers Numpy ###################################################
    mC = methodes_2020.MeshPlat(mC,direction)
    Namechamp = MEDchamp.getName()+direction+"Mean"
    MEDchampPlat = methodes_2020.NumpyToMEDField(champPlat,mC,temps,Namechamp)

    return (champPlat, MEDchampPlat)


def TimeAvg(MEDFileName,FieldName):

    field,mesh, TypeOfField = methodes_2020.MEDFileToMEDField_MultiTS(MEDFileName, FieldName)
    it = field.getIterations()

    # print(it)
    # k=field.getFieldOnMeshAtLevel(TypeOfField,it[0][0],it[0][1],0,mesh,0)
    # print(k.getNumberOfTuples())
    # print(mesh)
    print("phase 1 complete")


    ####################################################################
    ########################                    ########################
    ############      COEUR DE LA MOYENNE TEMPORELLE     ###############
    ########################                    ########################
    ####################################################################

    champMean = 0
    for i,ts in enumerate(it[1:]):
        ################################################################
        # REMARQUE : quand   ts  =  it[1] on a en fait i = 0 donc
        #            it[i]       =  it[0]
        #     ainsi (ts - it[i]) = (it[i] - it[i-1]) a tout instant
        
        dt = (ml.GetTimeAttachedOnFieldIteration(MEDFileName,FieldName,ts[0],ts[1])
              - ml.GetTimeAttachedOnFieldIteration(MEDFileName,FieldName,it[i][0],it[i][1])
              )
        
        champMean += ((field.getFieldOnMeshAtLevel(TypeOfField,ts[0],ts[1],0,mesh,0)
                    + field.getFieldOnMeshAtLevel(TypeOfField,it[i][0],it[i][1],0,mesh,0)
                    ) * dt/2 )

    ####################################################################
    ##################### RECUPERATION DU TEMPS ########################  
    # ATTN : Si le pas de temps n'est pas constant on pourra utiliser :
    # BIEN : MEDCoupling::GetTimeAttachedOnFieldIteration 	(
    # const std::string &  	fileName,
    # const std::string &  	fieldName,
    # int  	iteration,
    # int  	order 
    # )
    # NULL : MEDFileAnyTypeFieldMultiTS::getTimeStep 	( 	int  	iteration,	int  	order) 

    DT =  ml.GetTimeAttachedOnFieldIteration(MEDFileName,FieldName,it[-1][0],it[-1][1])
    DT -= ml.GetTimeAttachedOnFieldIteration(MEDFileName,FieldName,it[0][0],it[0][1])

    #########################                   ########################    
    ####################################################################

    
    ####################################################################
    ### ECRITURE DU NOUVEAU CHAMP 
    ####################################################################

    print("Time Span :", DT,"s")
    champMean = champMean/DT
    print("After mean:",champMean.getNumberOfTuples())
    
    MEDMean = ml.MEDCouplingFieldDouble(TypeOfField)
    MEDMean.setNature(ml.IntensiveMaximum)
    MEDMean.setName(FieldName+"_TAvg")
    MEDMean.setMesh(champMean.getMesh())
    MEDMean.setArray(champMean.getArray())
    ####################################################################
    
    ####################################################################
    ###### Conversion en Numpy !! 
    ####################################################################

    f2,mC = methodes_2020.en_tete(champMean)
    xn,yn,zn,champMean = methodes_2020.MEDFieldToNumpy(f2,mC)
    ####################################################################
    
    print(champMean.shape)

    return (champMean, MEDMean)


def Rij(MEDFileName,VelocityName):

    ####################################################################
    # Préparation
    FVelocity,FMesh,TOF = methodes_2020.MEDFileToMEDField_MultiTS(MEDFileName,VelocityName)
    Mesh = FMesh.getMeshAtLevel(0)
    it = FVelocity.getIterations()

        
    ####################################################################
    # Récupération des moyennes temporelles des vitesses
    # On dit ici que moyenne statistique égale moyenne temporelle (ergodicité)
    NPMeanVelocity, MEDMeanVelocity = outils_physique.TimeAvg(MEDFileName,VelocityName)

    ####################################################################
    # Preparation pour le MEDFile
    NewMEDFileName = MEDFileName.split('/')[-1]     # get only file name if path given
    NewMEDFileName = "Rij"+NewMEDFileName
    ml.WriteUMesh(NewMEDFileName,Mesh,True)  #SI ICI ON N'A PAS UN UMESH, REMPLACER PAR WriteMesh

    ####################################################################
    # Preparation pour les MEDchamps
        
    MEDFluctVel = ml.MEDCouplingFieldDouble.New(TOF,ml.ONE_TIME)
    MEDFluctVel.setName("Fluctuating Velocity")
    MEDFluctVel.setMesh(Mesh)
    
    MEDRij = ml.MEDCouplingFieldDouble.New(TOF,ml.ONE_TIME)
    MEDRij.setName("Reynolds stress")
    MEDRij.setMesh(Mesh)

    for i in it:
        t = ml.GetTimeAttachedOnFieldIteration(MEDFileName,VelocityName,i[0],i[1])
        temps = [t,i[0],i[1]]
        ####################################################################
        # Construction des vitesses fluctuantes
        # ATTN : Malgré son nom, NPVelocity est un MEDchamp à sa première ligne d'existence !
        #        Dès sa deuxième ligne d'existence, il est transformé en NumPy
        
        NPVelocity = FVelocity.getFieldOnMeshAtLevel(TOF,i[0],i[1],0,FMesh,0)
        NPVelocity,mC = methodes_2020.en_tete(NPVelocity)
        xn,yn,zn,NPVelocity = methodes_2020.MEDFieldToNumpy(NPVelocity,mC)

        NPFluctVel = NPVelocity-NPMeanVelocity
        sh    = NPFluctVel.shape    #sh[-1] doit valoir NumberOfComponents...
        ####################################################################
        # Construction des tensions de Reynolds
        shRij = sh+(sh[-1],)        #concatenation de tuples
        NPRij = np.zeros((shRij))
        for i in range(sh[-1]):
            for j in range(sh[-1]):
                NPRij[...,i,j] = NPFluctVel[...,i]*NPFluctVel[...,j]


        ####################################################################
        # Conversion en MEDchamp
        MEDFluctVel = methodes_2020.NumpyToMEDField_VRij(NPFluctVel,MEDFluctVel,mC,temps)
        MEDRij      = methodes_2020.NumpyToMEDField_VRij(NPRij     ,MEDRij,     mC,temps)
        
        ####################################################################
        # Ecriture dans le fichier MED des fluctuations et des tensions
        ml.WriteFieldUsingAlreadyWrittenMesh(NewMEDFileName,MEDFluctVel)
        ml.WriteFieldUsingAlreadyWrittenMesh(NewMEDFileName,MEDRij     )

    print("Rij fini")
    return (NPRij, MEDRij)

    



