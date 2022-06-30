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

################################################################################################
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

#########################################################################    
########### REPRI DU STAGE DE M2 ########################################    
def getBoutsBoiteI(LG,LD,mode='boite'):
    """
    NE BOUCLE PAS SUR LES BULLES
    Renvoi sour forme de liste de chaines de caracteres les suites de commandes
    pour extraire les boites autour des bulles, en prenant en ocmpte les cas des
    boites qui sortent du domaine.
    Entrées:
      - LG : liste de trois bool [tf,tf,tf] disant si la boite sort 'a gauche' pour la direction x, pour la direction y et pour la direction z
      - LD : liste de trois bool [tf,tf,tf] disant si la boite sort 'a droite' pour la direction x, pour la direction y et pour la direction z
    Sorties:
      - COMMANDES : liste de commandes pour extraire la boite autour de la bulle
    """
    COMMANDES = []
    ########################################################################
    # Organisons un peu, ce sera utile
    xg,yg,zg = LG[0],LG[1],LG[2]
    xd,yd,zd = LD[0],LD[1],LD[2]
    X,Y,Z    = xd+xg,yd+yg,zd+zg
    G,D      = xg+yg+zg, xd+yd+zd
    Ou_ca = [X, Y, Z]
    if 2 in Ou_ca : 
        print("la boite à bulle sort à gauche et \
a droite pour la direction %s. On ne traite pas ce cas..."%(Ou_ca.index(2)))
        COMMANDES = None
    
    
    
    if mode != 'correction':
        # print("xd, yd, zd")
        # print("xg, yg, zg")
        # print(xd, yd, zd)
        # print(xg, yg, zg)
        # print("G,D")
        # print(G,D)
        # print("X,Y,Z")
        # print(X,Y,Z)
        # print("Ou_ca")
        # print(Ou_ca)
        pass
    ########################################################################
    # Rangeons un peu
    # Debut de tranche, tranche complete, Fin de tranche
    tx = ['Idxm:',"Idxm:IdxM",":IdxM"]
    ty = ['Idym:',"Idym:IdyM",":IdyM"]
    tz = ['Idzm:',"Idzm:IdzM",":IdzM"]
    
    # Si on se sert de cette fonction pour la correction initiale du champ,
    # on doit faire un petit changement sur la decoupe des tranches
    if mode == 'correction':
        # tx = ['1:',":",":1"]
        # ty = ['1:',":",":1"]
        # tz = ['1:',":",":1"]
        tx[1] = ":"
        ty[1] = ":"
        tz[1] = ":"
    ########################################################################
    # petites fcts utiles pour aller chercher la bonne decoupe dans les
    # t[x,y,z] de juste au dessus
    ########################################################################
    # pour i :   0        1        2
    # f1 renvoie 0,1,1 ou 1,0,1 ou 1,1,0
    # f2 renvoie [2,2,2] - f1
    def f1(i):
        return(list((1-np.eye(3))[i]))
    def f2(i):
        return(list((1+np.eye(3))[i]))
    # pour i,j : 1,2      0,2      0,1
    # g1 renvoie 1,0,0 ou 0,1,0 ou 0,0,1
    # g2 renvoie 1,0,2 ou 0,2,1 ou 0,1,2
    # g3 renvoie 1,2,0 ou 2,1,0 ou 2,0,1
    # g4 renvoie 1,2,2 ou 2,1,2 ou 2,2,1
    def g1(i,j):
        return(list(np.eye(3)[abs(i+j-3)]))
    def g2(i,j):
        l = [1,1,1]
        l[j] = 2; l[i] = 0
        return(l)
    def g3(i,j):
        l = [1,1,1]
        l[j] = 0; l[i] = 2
        return(l)
    def g4(i,j):
        return(list(2*np.ones(3)-np.eye(3)[abs(i+j-3)]))
    ########################################################################
    
    ########################################################################
    ######## COEUR DE FONCTION                             #################
    ########################################################################
    
    if G+D == 0:
        # O| On ne sort pas. Nulle part
        Commande = "tous_les_bouts[0] = Numpy[{0},{1},{2}]".format(tx[1],ty[1],tz[1])
        
        COMMANDES = [Commande]
    
    ###################################
    if G+D == 1:
        # A| On sort d'un seul bout
        # drc dit par quelle direction on sort (par x, par y ou par z)
        drc = Ou_ca.index(1)
        # print("drc")
        # print(drc)
        # les s_i diront quelles decoupes aller chercher dans les t[x,y,z]
        s1,s2,s3 = int(f1(drc)[0]),int(f1(drc)[1]),int(f1(drc)[2])
        s4,s5,s6 = int(f2(drc)[0]),int(f2(drc)[1]),int(f2(drc)[2])
        # print("tx[s1],ty[s2],tz[s3]")
        # print("tx[s4],ty[s5],tz[s6]")
        # print(tx[s1],ty[s2],tz[s3])
        # print(tx[s4],ty[s5],tz[s6])
        Commande1 = "un_bout[0] = Numpy[{0},{1},{2}]".format(tx[s1],ty[s2],tz[s3])
        Commande2 = "un_bout[1] = Numpy[{0},{1},{2}]".format(tx[s4],ty[s5],tz[s6])
        
        Commande3 = "tous_les_bouts[0] = np.concatenate((un_bout[0],un_bout[1]),axis=%s)"%(drc)    
                    
                    
        COMMANDES = [Commande1,Commande2,Commande3]
    ####################################
    if G+D == 2:
        # B| On sort par deux cotes. 
        #    RAPPEL : Interdit de sortir a gauche ET a droite en meme temps
        #    --> les combinaisons possibles sont :
        #    drc1, drc2 = (1,2) ou (0,1) ou (0,2)
        drc1 = Ou_ca.index(1)
        drc2 = Ou_ca.index(1,drc1+1)
        # print("drc1,drc2")
        # print(drc1,drc2)
        # les s_i diront quelles decoupes aller chercher dans les t[x,y,z]
        s1, s2, s3  = int(g1(drc1,drc2)[0]),int(g1(drc1,drc2)[1]),int(g1(drc1,drc2)[2])
        s4, s5, s6  = int(g2(drc1,drc2)[0]),int(g2(drc1,drc2)[1]),int(g2(drc1,drc2)[2])
        s7, s8, s9  = int(g3(drc1,drc2)[0]),int(g3(drc1,drc2)[1]),int(g3(drc1,drc2)[2])
        s10,s11,s12 = int(g4(drc1,drc2)[0]),int(g4(drc1,drc2)[1]),int(g4(drc1,drc2)[2])
        
        Commande1 = "un_bout[0] = Numpy[{0},{1},{2}]".format(tx[s1],ty[s2],tz[s3])
        Commande2 = "un_bout[1] = Numpy[{0},{1},{2}]".format(tx[s4],ty[s5],tz[s6])
        Commande3 = "un_bout[2] = Numpy[{0},{1},{2}]".format(tx[s7],ty[s8],tz[s9])
        Commande4 = "un_bout[3] = Numpy[{0},{1},{2}]".format(tx[s10],ty[s11],tz[s12])
        
        Commande5 = "deux_bouts[0]   = np.concatenate((un_bout[0],un_bout[2]),axis=%s)"%(drc1) 
        Commande6 = "deux_bouts[1]   = np.concatenate((un_bout[1],un_bout[3]),axis=%s)"%(drc1) 
        Commande7 = "tous_les_bouts[0] = np.concatenate((deux_bouts[0],deux_bouts[1]),axis=%s)"%(drc2)
        
        COMMANDES = [Commande1,Commande2,Commande3,Commande4,Commande5,
                      Commande6,Commande7]
        
    ######################################
    if G+D ==3:
        # C| On sort par 3 cotes.
        #    RAPPEL : Interdit de sortir a gauche ET a droite en meme temps
        #    --> On sort donc forcement par x, par y et par z
        
        Commande1 = "un_bout[0] = Numpy[{0},{1},{2}]".format(tx[0],ty[0],tz[0])
        Commande2 = "un_bout[1] = Numpy[{0},{1},{2}]".format(tx[0],ty[0],tz[2])
        Commande3 = "un_bout[2] = Numpy[{0},{1},{2}]".format(tx[0],ty[2],tz[0])
        Commande4 = "un_bout[3] = Numpy[{0},{1},{2}]".format(tx[0],ty[2],tz[2])
        Commande5 = "un_bout[4] = Numpy[{0},{1},{2}]".format(tx[2],ty[0],tz[0])
        Commande6 = "un_bout[5] = Numpy[{0},{1},{2}]".format(tx[2],ty[0],tz[2])
        Commande7 = "un_bout[6] = Numpy[{0},{1},{2}]".format(tx[2],ty[2],tz[0])
        Commande8 = "un_bout[7] = Numpy[{0},{1},{2}]".format(tx[2],ty[2],tz[2])
    
        Commande9  = "deux_bouts[0] = np.concatenate((un_bout[0],un_bout[4]),axis=0)"
        Commande10 = "deux_bouts[1] = np.concatenate((un_bout[1],un_bout[5]),axis=0)"
        Commande11 = "deux_bouts[2] = np.concatenate((un_bout[2],un_bout[6]),axis=0)"
        Commande12 = "deux_bouts[3] = np.concatenate((un_bout[3],un_bout[7]),axis=0)"
    
        Commande13 = "quatre_bouts[0] = np.concatenate((deux_bouts[0],deux_bouts[2]),axis=1)"
        Commande14 = "quatre_bouts[1] = np.concatenate((deux_bouts[1],deux_bouts[3]),axis=1)"
    
        Commande15 = "tous_les_bouts[0] = np.concatenate((quatre_bouts[0],quatre_bouts[1]),axis=2)"
        
        COMMANDES = [Commande1 ,Commande2, Commande3, Commande4, Commande5,
                     Commande6 ,Commande7 ,Commande8 ,Commande9, Commande10,
                     Commande11,Commande12,Commande13,Commande14,Commande15]
        
    
    ##########################################
    
    ### ON NETTOIE LA MEMOIRE #########################################
    del(G,D,X,Y,Z,tx,ty,tz,Ou_ca,
        xg,yg,zg,
        xd,yd,zd)
    return(COMMANDES)
#########################################################################    
#########################################################################    
    


class Champ():
    """
    Classe de l'extraction d'un champ scalaire en 3D depuis le fichier curve de visit
    """
    def __init__(self):
        self.liste_nom_fichier="TEST_curve.curve"
        self.nb,self.data = 0,0 #np.loadtxt(self.nom_fichier, dtype=float, comments='#', delimiter=' ', unpack=True,
                             #    converters=None, skiprows=2, usecols=(0,1),  max_rows=None, like=None)
        self.L=0.004
        self.N=40 
        self.M=0 #int(self.N)
        self.field_kji=0#self.data.reshape(2*self.N,2*self.N,2*self.N)
        self.label_nom_du_champ="nom_du_champ_format_latex"
        self.lata_nom_du_champ="nom_du_champ_format_lata"
        self.temps=[]
        self.liste_data=[]
        self.liste_nb=[]
        self.liste_field_kji=[]
        self.liste_temps=[]
        
        
    def prepare_myself(self,liste_noms,jdd_name,repr_file):
        self.liste_nom_fichier=liste_noms
        self.liste_nom_fichier.sort(key=natural_keys)
        self.jdd = jdd_name
        self.repr_file = repr_file
        #self.nb,self.data = np.loadtxt(self.nom_fichier, dtype=float, comments='#', delimiter=' ', unpack=True,
                                       #converters=None, skiprows=2, usecols=(0,1),  max_rows=None, like=None)
        # self.L=0.004
        self.L=dtool.getParam(self.jdd+".data","uniform_domain_size_i")
        # self.N=40 
        self.N=int(dtool.getParam(self.jdd+".data","nbelem_i"))
        self.M=int(self.N)
        #self.field_kji=self.data.reshape(2*self.N,2*self.N,2*self.N)
        self.label_nom_du_champ="nom_du_champ_format_latex"                               
        self.lata_nom_du_champ="nom_du_champ_format_lata"
        self.liste_temps = []
        for nom_fichier in self.liste_nom_fichier:
            f=open(nom_fichier)
            time = float(f.readline().split('TIME')[-1])
            self.liste_temps.append(time)
                
    def translate_noms(self,read_name_from_file=True):
        """
        Lis le nom du champ au format lata si read_name_from_file=True
        Traduit ce nom lata en nom bien ecrit pour les label matplotlib
        """
        if read_name_from_file:
            f = open(self.nom_fichier,"r")
            l = f.readline(); l = f.readline()
            f.close()
            self.lata_nom_du_champ = l.lstrip('# ').rstrip('\n').rstrip('_dual').rstrip('_DOM').lower()
        self.label_nom_du_champ = self.lata_nom_du_champ.rstrip('\n').rstrip('_dual').rstrip('_DOM').rstrip('_ELEM').rstrip('_FACES').lower().split('_')
        self.label_nom_du_champ = r"$"+self.label_nom_du_champ[0]+"_"+self.label_nom_du_champ[-1]+"$"
        
    ####################################################################        
    def load_fields(self):
        
        # self.liste_data,self.liste_nb,self.liste_field_kji,self.liste_temps = [],[],[],[]
        for nom_fichier in self.liste_nom_fichier:
            
            nb,data = np.loadtxt(nom_fichier, dtype=float, comments='#', delimiter=' ', unpack=True,
                                       converters=None, skiprows=2, usecols=(0,1),  max_rows=None, like=None)
            if 'dual' in nom_fichier:
                field_kji = data.reshape(2*self.N,2*self.N,2*self.N)
            else:
                field_kji = data.reshape(self.N,self.N,self.N)
            
            self.liste_data.append(data)
            self.liste_nb.append(nb)
            self.liste_field_kji.append(field_kji)
    
    ####################################################################        
    def load_fields_with_raw(self):
        
        # self.liste_data,self.liste_nb,self.liste_field_kji,self.liste_temps = [],[],[],[]
        for nom_fichier in self.liste_nom_fichier:
            
            nb,data = np.fromfile(nom_fichier, dtype=float, comments='#', delimiter=' ', unpack=True,
                                       converters=None, skiprows=2, usecols=(0,1),  max_rows=None, like=None)
            if 'dual' in nom_fichier:
                field_kji = data.reshape(2*self.N,2*self.N,2*self.N)
            else:
                field_kji = data.reshape(self.N,self.N,self.N)
                
            self.liste_data.append(data)
            self.liste_nb.append(nb)
            self.liste_field_kji.append(field_kji)
    ####################################################################        
    def elem_field_to_face_field(self):
        """
        Interpole un champ aux elements en un champ aux faces
         - La premiere face a gauche est interieure au domaine
         - La derniere face a droite est exterieur, construite par interpolation
         - Retourne une liste 
        """
        liste_x_faces_field_kji,liste_y_faces_field_kji,liste_z_faces_field_kji = self.liste_field_kji.copy(),self.liste_field_kji.copy(),self.liste_field_kji.copy()
        for it,field_kji in enumerate(self.liste_field_kji):
            # Faces interieures au domaine. La prmeiere face correspond à "i+1/2".
            liste_x_faces_field_kji[it][...,:,:,:-1] = field_kji[...,:,:,1:]+field_kji[...,:,:,:-1]
            liste_y_faces_field_kji[it][...,:,:-1,:] = field_kji[...,:,1:,:]+field_kji[...,:,:-1,:]
            liste_z_faces_field_kji[it][...,:-1,:,:] = field_kji[...,1:,:,:]+field_kji[...,:-1,:,:]
            # Derniere face "a droite", periodicite du domaine
            liste_x_faces_field_kji[it][...,:,:,-1] = field_kji[...,:,:,0] +field_kji[...,:,:,-1]
            liste_y_faces_field_kji[it][...,:,-1,:] = field_kji[...,:,0,:] +field_kji[...,:,-1,:]
            liste_z_faces_field_kji[it][...,-1,:,:] = field_kji[...,0,:,:] +field_kji[...,-1,:,:]
            return(liste_x_faces_field_kji,liste_y_faces_field_kji,liste_z_faces_field_kji)
            
    def get_bubble_centers(self):
        diphasique = True
        if diphasique:
            #_,_,_,_,self.alpha,self.nb_bulles,_,_,_,_,_,_ = dtool.get_prop(jdd=self.jdd+".data", repr_file=self.repr_file, Ret=0.)
            cbx,cby,cbz = glob.glob("OUT/"+self.jdd+"_bulles_centre_x.out"),glob.glob("OUT/"+self.jdd+"_bulles_centre_y.out"),glob.glob("OUT/"+self.jdd+"_bulles_centre_z.out")
            centres_x = np.loadtxt(cbx[0], dtype=float)
            centres_y = np.loadtxt(cby[0], dtype=float)
            centres_z = np.loadtxt(cbz[0], dtype=float)
            
            coords = np.array([centres_x[:,1:],centres_y[:,1:],centres_z[:,1:]])
            self.liste_temps_vitesse_bulle, self.liste_vitesses_bulles = dtool.getVelocity(coords,centres_x[:,0],centres_x.shape[-1],self.jdd+".data")
            
            self.liste_centres_x = np.zeros((len(self.liste_temps),len(centres_x[0,1:])))
            self.liste_centres_y = np.zeros((len(self.liste_temps),len(centres_y[0,1:])))
            self.liste_centres_z = np.zeros((len(self.liste_temps),len(centres_z[0,1:])))
            # self.liste_vitesses_bulles = np.zeros((len(self.liste_temps)))
            for it, t in enumerate(self.liste_temps):
                good_line = np.argmin(abs(centres_x[:,0]-t))
                self.liste_centres_x[it,:] = centres_x[good_line,1:]
                self.liste_centres_y[it,:] = centres_y[good_line,1:]
                self.liste_centres_z[it,:] = centres_z[good_line,1:]
                # self.liste_vitesses_bulles[it] = liste_vitesses_bulles[good_line]
            # le numpy self.liste_centres_xyz est agence comme : [iteration,numero_bulle,kji]   
            self.liste_centres_xyz = np.array([self.liste_centres_x,self.liste_centres_y,self.liste_centres_z]).transpose(1,2,0)
            
            del(centres_x,centres_y,centres_z)
        else:
            print("Pas de bulles : n'appelle pas get_bubble_centers !")            
            
            
    def get_bubble_boxes(self):  
        diphasique = True    
        if diphasique:    
            # Parametres physique des bulles
            self.vb = dtool.getParam(self.jdd+".data", "vol_bulle_monodisperse")
            self.db = 2.*pow(3./(4*mt.pi)*self.vb,1./3.)
            
            # Extentsions arbitraires autour du centre des bulles
            trailing = 1.*self.db
            front = 1.*self.db
            left_hand = 1.*self.db 
            
            # Conversion en nombre de cases (en indices)
            self.left_hand_in_index, self.left_hand_in_index, self.front_in_index = self.pos_to_kji((front, left_hand, left_hand))
            self.left_hand_in_index, self.left_hand_in_index, self.trailing_in_index = self.pos_to_kji((trailing, left_hand, left_hand))
            if 'dual' in self.liste_nom_fichier[0]:
                # Retranche la moitie du nombre de cellules pour mettre le 0,0,0 au milieu
                self.front_in_index, self.left_hand_in_index, self.left_hand_in_index = self.front_in_index-self.N, self.left_hand_in_index-self.N, self.left_hand_in_index-self.N
                self.trailing_in_index = self.trailing_in_index-self.N
            else:
                # Retranche la moitie du nombre de cellules pour mettre le 0,0,0 au milieu
                self.front_in_index, self.left_hand_in_index, self.left_hand_in_index = int(self.front_in_index-self.N/2), int(self.left_hand_in_index-self.N/2), int(self.left_hand_in_index-self.N/2)
                self.trailing_in_index = int(self.trailing_in_index-self.N/2)
            
            
            # Conversion de la position du centre des bulles xyz en indice kji       [temps,bulle,kji]
            self.liste_centres_kji = self.liste_time_bubble_pos_to_liste_time_bubble_kji(self.liste_centres_xyz)
            self.liste_bubble_field_kji = np.zeros((self.liste_centres_kji.shape[0],
                                                    self.liste_centres_kji.shape[1],
                                                    2*self.left_hand_in_index,
                                                    2*self.left_hand_in_index,
                                                    self.front_in_index+self.trailing_in_index))
            # Obtention des indices des bornes des boites                            [temps,bulle,kji]
            self.liste_minima_boites_kji = np.array([self.liste_centres_kji[...,0]-self.left_hand_in_index,
                                                     self.liste_centres_kji[...,1]-self.left_hand_in_index,
                                                     self.liste_centres_kji[...,2]-self.trailing_in_index]).transpose(1,2,0)
            self.liste_Maxima_boites_kji = np.array([self.liste_centres_kji[...,0]+self.left_hand_in_index,
                                                     self.liste_centres_kji[...,1]+self.left_hand_in_index,
                                                     self.liste_centres_kji[...,2]+self.front_in_index]).transpose(1,2,0)
            
            
            # Reponds a la question "Combien de fois sort-on a gauche/ a droite", avec "a gauche" qui signifie "par valeur inferieure", implicitement          
            # On met le 2*N si on travail sur un champ du domaine dual. Autrement on prend juste N !!!
            cumul_Droite_Sort_kji,cumul_Gauche_Sort_kji = 0*self.liste_Maxima_boites_kji,0*self.liste_minima_boites_kji
            self.liste_corrigee_minima_boites_kji,self.liste_corrigee_Maxima_boites_kji = self.liste_minima_boites_kji,self.liste_Maxima_boites_kji
            if 'dual' in self.liste_nom_fichier[0]:
                while (self.liste_corrigee_minima_boites_kji<0).any() :
                    cumul_Gauche_Sort_kji += 1*self.liste_corrigee_minima_boites_kji<0
                    self.liste_corrigee_minima_boites_kji = self.liste_corrigee_minima_boites_kji + 2*self.N*(self.liste_corrigee_minima_boites_kji<0)
                while (self.liste_corrigee_Maxima_boites_kji>2*self.N).any() :
                    cumul_Droite_Sort_kji += 1*self.liste_corrigee_Maxima_boites_kji>2*self.N 
                    self.liste_corrigee_Maxima_boites_kji = self.liste_corrigee_Maxima_boites_kji - 2*self.N*(self.liste_corrigee_Maxima_boites_kji>2*self.N)

            else:
                while (self.liste_corrigee_minima_boites_kji<0).any() :
                    cumul_Gauche_Sort_kji += 1*self.liste_corrigee_minima_boites_kji<0
                    self.liste_corrigee_minima_boites_kji = self.liste_corrigee_minima_boites_kji + self.N*(self.liste_corrigee_minima_boites_kji<0)
                while (self.liste_corrigee_Maxima_boites_kji>self.N).any() :
                    cumul_Droite_Sort_kji += 1*self.liste_corrigee_Maxima_boites_kji>self.N 
                    self.liste_corrigee_Maxima_boites_kji = self.liste_corrigee_Maxima_boites_kji - self.N*(self.liste_corrigee_Maxima_boites_kji>self.N)
            
            
            # Reponds a la question "Sort-on a gauche/ a droite", avec "a gauche" qui signifie "par valeur inferieure", implicitement          
            # On met le 2*N si on travail sur un champ du domaine dual. Autrement on prend juste N !!!
            if 'dual' in self.liste_nom_fichier[0]:
                Gauche_Sort_kji,Droite_Sort_kji = self.liste_minima_boites_kji<0 , self.liste_Maxima_boites_kji>2*self.N 
            else:
                Gauche_Sort_kji,Droite_Sort_kji = self.liste_minima_boites_kji<0 , self.liste_Maxima_boites_kji>self.N 
            # commandes_extraction : [it][ib,1,longueur_variable], ou longueur_variable vaut :
            #  - 1  : la boite ne sort pas du domaine
            #  - 3  : la boite sort pour une des directions
            #  - 7  : la boite sort pour deux des troix directions
            #  - 15 : la boite sort dans toutes les directions
            self.liste_commandes_extraction = self.getBoutsBoite_liste(Gauche_Sort_kji,Droite_Sort_kji,cumul_Gauche_Sort_kji,cumul_Droite_Sort_kji)
            
                            
        else:
            print("Pas de bulles : n'appelle pas get_field_around_bubble !")            
       
    
    def getBoutsBoite_liste(self,LG_kji,LD_kji,cumul_LG_kji,cumul_LD_kji,mode='boite'):
        """
        Renvoi sour forme de liste de chaines de caracteres les suites de commandes
        pour extraire les boites autour des bulles, en prenant en ocmpte les cas des
        boites qui sortent du domaine.
        Entrées:
        - LG_kji : array de bool de taille [temps,bulle,3] disant si la boite de la bulle, a l'instant temps, sort 'a gauche' pour la direction z, pour la direction y et pour la direction x
        - LD_kji : array de bool de taille [temps,bulle,3] disant si la boite de la bulle, a l'instant temps, sort 'a droite' pour la direction z, pour la direction y et pour la direction x
        Sorties:
        - COMMANDES : liste de commandes pour extraire la boite autour de la bulle
        
        REMARQUES : Deux nouveaux cas de figures foncitonnent pour le setup numerique habituel. 
                    Il faut les generaliser pour que ces deux cas foncitonnent quel que soit le setup (les parties soulignees sont celles a generaliser)
                      - La boite sort a gauche ET a droite dans la direciton x uniquement
                                                                --------------
                      - La boite sort de plus d'une longueur de boite par la gauche dans la direction x uniquement 
                                                                      -------------      --------------
        """
        ## Preparation des commandes
        # liste pour stocker les commander d'extraction
        COMMANDES = []
        # nom du champ. Dans cette casse il s'appelle 'self.liste_field_kji', puis on appelera la commande à chaque pas de temps 
        nom_du_champ_format_numpy = 'self.liste_field_kji[it]'
        # nom des indices qui referent aux min/Max des boites qu'on veux extraire. Dans cette classe, il s'appallent liste_minima_boites_kji/liste_Maxima_boites_kji, 
        #   puis on appelera la commande à chaque pas de temps et pour chaque bulle
        nom_d_indices_mini_boites_kji = 'self.liste_corrigee_minima_boites_kji[it,ib]'
        nom_d_indices_Maxi_boites_kji = 'self.liste_corrigee_Maxima_boites_kji[it,ib]'
        ########################################################################
        # Organisons un peu, ce sera utile
        xg,yg,zg = 1*LG_kji[...,2],1*LG_kji[...,1],1*LG_kji[...,0]
        xd,yd,zd = 1*LD_kji[...,2],1*LD_kji[...,1],1*LD_kji[...,0]
        X,Y,Z    = xd+xg,yd+yg,zd+zg
        G,D      = 1*(xg+yg+zg), 1*(xd+yd+zd)
        # HACK : " 1* " -> pour transformer un bool en int
        GD = 1*np.sum(LG_kji,axis=-1,keepdims=True)+1*np.sum(LD_kji,axis=-1,keepdims=True)
        Ou_ca = np.array([X, Y, Z],dtype=int).transpose(1,2,0)
        # print("xg,yg,zg",xg,yg,zg)
        # print("xd,yd,zd",xd,yd,zd)
        nb_bulle, nb_temps = Ou_ca.shape[1], Ou_ca.shape[0]
        # Objet qui va contenir les suites d'instructions pour extraire la fameuse
        COMMANDES = np.array([[[None]] * nb_bulle] * nb_temps)
        
        ########################################################################
        # Rangeons un peu
        # ti = [ Debut de tranche, tranche complete, Fin de tranche ]
        ########################################################################
        # tx,ty,tz sont ok mais cest dans les bout qu'on met pas comme il faut !!!
        ########################################################################
        ti = [
              'int(np.round('+str(nom_d_indices_mini_boites_kji.replace('ib','ib,2'))+',0))'+':',
              'int(np.round('+str(nom_d_indices_mini_boites_kji.replace('ib','ib,2'))+',0))'+
          ":"+'int(np.round('+str(nom_d_indices_Maxi_boites_kji.replace('ib','ib,2'))+',0))',
          ":"+'int(np.round('+str(nom_d_indices_Maxi_boites_kji.replace('ib','ib,2'))+',0))'
             ]                                                                       
        tj = [                                                                       
              'int(np.round('+str(nom_d_indices_mini_boites_kji.replace('ib','ib,1'))+',0))'+':',
              'int(np.round('+str(nom_d_indices_mini_boites_kji.replace('ib','ib,1'))+',0))'+
          ":"+'int(np.round('+str(nom_d_indices_Maxi_boites_kji.replace('ib','ib,1'))+',0))',
          ":"+'int(np.round('+str(nom_d_indices_Maxi_boites_kji.replace('ib','ib,1'))+',0))'
             ]                     
        tk = [                     
              'int(np.round('+str(nom_d_indices_mini_boites_kji.replace('ib','ib,0'))+',0))'+':',
              'int(np.round('+str(nom_d_indices_mini_boites_kji.replace('ib','ib,0'))+',0))'+
          ":"+'int(np.round('+str(nom_d_indices_Maxi_boites_kji.replace('ib','ib,0'))+',0))',
          ":"+'int(np.round('+str(nom_d_indices_Maxi_boites_kji.replace('ib','ib,0'))+',0))'
             ]
        ########################################################################
        # tx,ty,tz sont ok mais cest dans les bout qu'on met pas comme il faut !!!
        ########################################################################
             
        # Si on se sert de cette fonction pour la correction initiale du champ,
        # on doit faire un petit changement sur la decoupe des tranches
        if mode == 'correction':
            ti[1] = ":"
            tj[1] = ":"
            tk[1] = ":"
        ########################################################################
        # petites fcts utiles pour aller chercher la bonne decoupe dans les
        # t[x,y,z] de juste au dessus
        ########################################################################
        # pour i :   0        1        2
        def f1(i):
            """ f1 renvoie 0,1,1 ou 1,0,1 ou 1,1,0 """
            return(list((1-np.eye(3))[i]))
        def f2(i):
            """ f2 renvoie [2,2,2] - f1 """
            return(list((1+np.eye(3))[i]))
        # pour i,j : 1,2      0,2      0,1
        def g1(i,j):
            """ g1 renvoie 1,0,0 ou 0,1,0 ou 0,0,1 """
            return(list(np.eye(3)[abs(i+j-3)]))
        def g2(i,j):
            """ g2 renvoie 1,0,2 ou 0,2,1 ou 0,1,2 """
            l = [1,1,1]
            l[j] = 2; l[i] = 0
            return(l)
        def g3(i,j):
            """ g3 renvoie 1,2,0 ou 2,1,0 ou 2,0,1 """
            l = [1,1,1]
            l[j] = 0; l[i] = 2
            return(l)
        def g4(i,j):
            """ g4 renvoie 1,2,2 ou 2,1,2 ou 2,2,1 """
            return(list(2*np.ones(3)-np.eye(3)[abs(i+j-3)]))
        ########################################################################
        
        ########################################################################
        ######## COEUR DE FONCTION                             #################
        ########################################################################
        
        for it in range(nb_temps):
            for ib in range(nb_bulle):
                if GD[it,ib] == 0:
                    # O| On ne sort pas. Nulle part.
                    # Commande = "tous_les_bouts[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tx[1],ty[1],tz[1])
                    Commande = "tous_les_bouts[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[1],tj[1],ti[1])
                    COMMANDES[it,ib,0] = [Commande]
                    
                elif GD[it,ib] == 1:
                    # A| On sort d'un seul bout
                    # drc dit par quelle direction on sort (par x, par y ou par z)
                    drc = np.argwhere(Ou_ca[it,ib,:]==1)[0][0]
                    # les s_i diront quelles decoupes aller chercher dans les t[x,y,z]
                    s1,s2,s3 = int(f1(drc)[0]),int(f1(drc)[1]),int(f1(drc)[2])
                    s4,s5,s6 = int(f2(drc)[0]),int(f2(drc)[1]),int(f2(drc)[2])
                    
                    if not( (cumul_LG_kji[it,ib,:]>1).any() or (cumul_LD_kji[it,ib,:]>1).any() ):
                        # Si on sort de moins d'un domaine complet :
                        Commande1 = "un_bout[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s3],tj[s2],ti[s1])
                        # Commande1 = "un_bout[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tx[s1],ty[s2],tz[s3])
                        Commande2 = "un_bout[1] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s6],tj[s5],ti[s4])
                        # Commande2 = "un_bout[1] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tx[s4],ty[s5],tz[s6])
                        
                        # Commande3 = "tous_les_bouts[0] = np.concatenate((un_bout[0],un_bout[1]),axis=%s)"%(drc)    
                        # drc=0 correspond a la direction x, ce qui correspond a axis=2 en convention kji
                        Commande3 = "tous_les_bouts[0] = np.concatenate((un_bout[0],un_bout[1]),axis=%s)"%(2-drc)   
                                                            
                        COMMANDES[it,ib,0] = [Commande1,Commande2,Commande3]
                    else:
                        # Si on sort de plus d'un domaine complet FONCTIONNE POUR LA DIRECTION X, en sortant a gauche UNIQUEMENT
                        Commande1 = "un_bout[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s3],tj[s2],ti[s1])
                        Commande2 = "un_bout[1] = np.concatenate(([{0}[1:-1,1:-1,:]]*{1}),axis={2})".format(nom_du_champ_format_numpy,int(cumul_LG_kji[it,ib,2]-1),"2")
                        Commande3 = "un_bout[2] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s6],tj[s5],ti[s4])
                        
                        Commande4 = "tous_les_bouts[0] = np.concatenate((un_bout[0],un_bout[1],un_bout[2]),axis=%s)"%(2-drc)  
                        COMMANDES[it,ib,0] = [Commande1,Commande2,Commande3,Commande4]
                
                elif GD[it,ib] == 2:
                    # B| On sort par deux cotes. 
                    #    RAPPEL : Interdit de sortir a gauche ET a droite en meme temps SAUF pour la direction x (direction de la gravite)
                    #    --> les combinaisons possibles sont :
                    #    drc1, drc2 = (1,2) ou (0,1) ou (0,2) si on ne sort pas des deux bouts par x
                    #    drc1, drc2 = (0,0) si on sort des deux bouts par x
                    ########################################################
                    # Foncitonne entierement que pour g selon x / que si on sort des 2 cotes selon x ...
                    # A ADAPTER POUR QUE TOUT FONCITONNE AUSSI BIEN SELON
                    # TOUTES LES DIRECTION
                    ########################################################
                    if Ou_ca[it,ib,0] == 2:
                        # Les deux sorties se font selon la direction x
                        # Dans les deux autres directions, on ne sort pas
                        drc1,drc2 = 0,0
                        # s1, s2, s3  = int(g1(drc1,drc2)[0]),int(g1(drc1,drc2)[1]),int(g1(drc1,drc2)[2])
                        # s4, s5, s6  = int(g2(drc1,drc2)[0]),int(g2(drc1,drc2)[1]),int(g2(drc1,drc2)[2])
                        # s7, s8, s9  = int(g3(drc1,drc2)[0]),int(g3(drc1,drc2)[1]),int(g3(drc1,drc2)[2])
                        
                        Commande1 = "un_bout[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[1],tj[1],ti[0])
                        Commande2 = "un_bout[1] = {0}[{1},{2},:]".format(nom_du_champ_format_numpy,tk[1],tj[1])
                        Commande3 = "un_bout[2] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[1],tj[1],ti[2])
                        
                        Commande7 = "tous_les_bouts[0] = np.concatenate((un_bout[0],un_bout[1],un_bout[2]),axis=%s)"%("2") 
                        
                        COMMANDES[it,ib,0] = [Commande1,Commande2,Commande3,Commande7]
                    else:
                        drc1,drc2 = np.argwhere(Ou_ca[it,ib,:]==1)[:,0]
                        # les s_i diront quelles decoupes aller chercher dans les t[x,y,z]
                        s1, s2, s3  = int(g1(drc1,drc2)[0]),int(g1(drc1,drc2)[1]),int(g1(drc1,drc2)[2])
                        s4, s5, s6  = int(g2(drc1,drc2)[0]),int(g2(drc1,drc2)[1]),int(g2(drc1,drc2)[2])
                        s7, s8, s9  = int(g3(drc1,drc2)[0]),int(g3(drc1,drc2)[1]),int(g3(drc1,drc2)[2])
                        s10,s11,s12 = int(g4(drc1,drc2)[0]),int(g4(drc1,drc2)[1]),int(g4(drc1,drc2)[2])
                        
                        Commande1 = "un_bout[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s3],tj[s2],ti[s1])
                        Commande2 = "un_bout[1] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s6],tj[s5],ti[s4])
                        Commande3 = "un_bout[2] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s9],tj[s8],ti[s7])
                        Commande4 = "un_bout[3] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[s12],tj[s11],ti[s10])
                        
                        Commande5 = "deux_bouts[0]   = np.concatenate((un_bout[0],un_bout[2]),axis=%s)"%(2-drc1) 
                        Commande6 = "deux_bouts[1]   = np.concatenate((un_bout[1],un_bout[3]),axis=%s)"%(2-drc1) 
                        Commande7 = "tous_les_bouts[0] = np.concatenate((deux_bouts[0],deux_bouts[1]),axis=%s)"%(2-drc2)
                        
                        COMMANDES[it,ib,0] = [Commande1,Commande2,Commande3,Commande4,Commande5,
                                                Commande6,Commande7]
                    
                elif GD[it,ib] == 3:
                    # C| On sort par 3 cotes.
                    #    RAPPEL : Interdit de sortir a gauche ET a droite en meme temps
                    #    --> On sort donc forcement par x, par y et par z
                    Commande1 = "un_bout[0] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[0],tj[0],ti[0])
                    Commande2 = "un_bout[1] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[0],tj[0],ti[2])
                    Commande3 = "un_bout[2] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[0],tj[2],ti[0])
                    Commande4 = "un_bout[3] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[0],tj[2],ti[2])
                    Commande5 = "un_bout[4] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[2],tj[0],ti[0])
                    Commande6 = "un_bout[5] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[2],tj[0],ti[2])
                    Commande7 = "un_bout[6] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[2],tj[2],ti[0])
                    Commande8 = "un_bout[7] = {0}[{1},{2},{3}]".format(nom_du_champ_format_numpy,tk[2],tj[2],ti[2])
                    
                    Commande9  = "deux_bouts[0] = np.concatenate((un_bout[0],un_bout[4]),axis=0)"
                    Commande10 = "deux_bouts[1] = np.concatenate((un_bout[1],un_bout[5]),axis=0)"
                    Commande11 = "deux_bouts[2] = np.concatenate((un_bout[2],un_bout[6]),axis=0)"
                    Commande12 = "deux_bouts[3] = np.concatenate((un_bout[3],un_bout[7]),axis=0)"
                    
                    Commande13 = "quatre_bouts[0] = np.concatenate((deux_bouts[0],deux_bouts[2]),axis=1)"
                    Commande14 = "quatre_bouts[1] = np.concatenate((deux_bouts[1],deux_bouts[3]),axis=1)"
                    
                    Commande15 = "tous_les_bouts[0] = np.concatenate((quatre_bouts[0],quatre_bouts[1]),axis=2)"
                    
                    COMMANDES[it,ib,0] = [Commande1 ,Commande2, Commande3, Commande4, Commande5,
                                            Commande6 ,Commande7 ,Commande8 ,Commande9, Commande10,
                                            Commande11,Commande12,Commande13,Commande14,Commande15]
                    pass
                pass
            pass
            
        return(COMMANDES) 
    
    
    
    def get_field_around_bubble(self):
        """
        A FINIR
        Objet : liste_bubble_field_kji ses indices sont : [t,b,k,j,i]
                     -> t : temps
                     -> b : bulle
                     -> k : dimension z
                     -> j : dimension y
                     -> i : dimension x
        """
        diphasique = True
        if diphasique:
            for it,field_kji in enumerate(self.liste_field_kji) : 
                commandes = self.liste_commandes_extraction[it]
                
                un_bout = [None]*8
                deux_bouts = [None]*4
                quatre_bouts = [None]*2
                tous_les_bouts = [None]
                
                
                
                for ib in range(self.liste_Maxima_boites_kji.shape[1]):
                    print(self.liste_bubble_field_kji[it,ib,...].shape)
                    for c in commandes[ib,0] : 
                        exec(c)
                    self.liste_bubble_field_kji[it,ib,...] = tous_les_bouts[0].copy()#.transpose(2,1,0).copy()
        else:
            print("Pas de bulles : n'appelle pas get_field_around_bubble !")
            
    def moyenne_temporelle(self,champ):
        """
        le champ doit etre au format [t,b,k,j,i] et rien d'autre
        """
        return(np.mean(champ,axis=0))
        
        
    ####################################################################
    def nbc_to_kji(self,nb):
        """
        Numero de cellule nb vers indice (i,j,k)
        D'abord i, ensuite j, enfin k.
        k -> Quotient D.E. de N par 80**2. Reste note R_ji
        j -> Quotient D.E. de R_ji par 80**1. Reste note R_i
        i -> R_i
        """
        if 'dual' in self.liste_nom_fichier[0]:
            k = int(nb/(2*self.N)**2)
            j = int((nb-k*(2*self.N)**2)/(2*self.N))
            i = int(nb-k*(2*self.N)**2 - j*(2*self.N))
        else:
            k = int(nb/(self.N)**2)
            j = int((nb-k*(self.N)**2)/(self.N))
            i = int(nb-k*(self.N)**2 - j*(self.N))
            
        return(k,j,i)
    
    def kji_to_nbc(self,kji):
        """ indice (i,j,k) vers numero de cellule nb"""
        if 'dual' in self.liste_nom_fichier[0]:
            return(kji[0] + kji[1]*(2*slef.N) + kji[2]*(2*slef.N)**2)
        else:
            return(kji[0] + kji[1]*(slef.N) + kji[2]*(slef.N)**2)
    
    def kji_to_pos(self,kji):
        """indice (i,j,k) vers position (x,y,z)"""
        i,j,k = kji[2],kji[1],kji[0]
        if 'dual' in self.liste_nom_fichier[0]:
            x = self.L/(2*self.N) * i - self.L/2
            y = self.L/(2*self.N) * j - self.L/2
            z = self.L/(2*self.N) * k - self.L/2
        else:
            i,j,k = kji[2],kji[1],kji[0]
            x = self.L/(self.N) * i - self.L/2
            y = self.L/(self.N) * j - self.L/2
            z = self.L/(self.N) * k - self.L/2
            
        return(x,y,z)
    
    def pos_to_kji(self,xyz):
        """position (x,y,z) vers indice (i,j,k)"""
        x,y,z = xyz[0], xyz[1], xyz[2]
        if 'dual' in self.liste_nom_fichier[0]:
            i = int((x+self.L/2)*2*self.N/self.L) 
            j = int((y+self.L/2)*2*self.N/self.L) 
            k = int((z+self.L/2)*2*self.N/self.L)
        else:
            # L ERREUR ELLE VIENT DE LA !!! 
            # REFAIRE LE PTIT DESSIN DE LA FCT AFFINE POUR VOIR QUEL FACTEUR
            # EST MAUVAIS !!!
            i = int((x+self.L/2)*self.N/self.L)             
            j = int((y+self.L/2)*self.N/self.L)             
            k = int((z+self.L/2)*self.N/self.L)            
            
        return(k,j,i)
        
    def liste_time_bubble_pos_to_liste_time_bubble_kji(self,xyz):
        """liste positions (iteration,numero_bulle,x,y,z) vers indice (iteration,numero_bulle,i,j,k)"""
        x,y,z = xyz[...,0], xyz[...,1], xyz[...,2]
        if 'dual' in self.liste_nom_fichier[0]:
            i = (x+self.L/2)*2*self.N/self.L 
            j = (y+self.L/2)*2*self.N/self.L 
            k = (z+self.L/2)*2*self.N/self.L
        else:
            i = (x+self.L/2)*self.N/self.L
            j = (y+self.L/2)*self.N/self.L
            k = (z+self.L/2)*self.N/self.L
            
        # transpose permet d'etre arrange comme : temps,bulle,indice
        return(np.array([k,j,i]).transpose(1,2,0))
                
    
    def nbc_to_pos(self,nb):
        """
        Nombre nb vers position (x,y,z)
        """
        return(self.kji_to_pos((self.nbc_to_kji(nb))))
    
    def pos_to_sph(self,xyz):
        """
        position cartesienne (x,y,z) vers postion spherique (r,theta,phi)
        /!\ La gravite est selon x, donc on choisi
            Theta : angle x,M
            Phi   : angle y,z
        """
        x,y,z = xyz[0], xyz[1], xyz[2]
        
        r = mt.sqrt(x*x + y*y + z*z)
        theta = mt.atan( mt.sqrt(y*y + z*z) / (x+1e-20))
        phi = mt.atan(z/(y+1e-20))
        
        return(r,theta,phi)
    
    def sph_to_pos(self,rtp):
        """
        Position spherique vers position cartesienne
        /!\ La gravite est selon x, donc on choisi
            Theta : angle x,M
            Phi   : angle y,z
        """
        r, theta, phi = rtp[0], rtp[1], rtp[2]
        x = cos(theta)
        y = sin(theta)*cos(phi)
        z = sin(theta)*sin(phi)
        
        return(r,theta,phi)
    
    ####################################################################
    def sph_mean(self):
        """
        Moyenne 3D du champ scalaire data, coquille spherique
        data(nbc) -> moy_{theta,phi}[data](r)
        """
        self.M=int(self.N)
        self.fieldR = np.zeros(self.M+1)
        self.R_max = (self.pos_to_sph(self.nbc_to_pos(0))[0]) / mt.sqrt(3) # pos_to_sph(nb_to_pos(0)) : cellule dans un COIN. On veut sur la frontiere, mais atteingnable pour tout x,y,z 
        self.dr = self.R_max/(self.M)
        self.volume_cellule = np.abs(self.nbc_to_pos(0)[0] - self.nbc_to_pos(1)[0]) **3 
        self.R = np.linspace(0,self.R_max,self.M+1)
        for i in self.nb:
            r = self.pos_to_sph(self.nbc_to_pos(i))[0]
            if r+self.dr>r and r>=r-self.dr and r<self.R_max:
                difference = np.abs(self.R - r) 
                indice_r = np.argwhere(difference == np.min(difference))[0][0]
                self.fieldR[indice_r] += self.data[int(i)] * self.volume_cellule
        
        self.fieldR/=(8*np.pi*(self.R+self.dr/2.)**2 *self.dr)
        # return(self.R,self.fieldR)
        
    def data_nbc_to_sph(self,rosette_board=False):
        """
        A partir du champ au format 1D, ressort le champ en spherique : 
        data(nbc) -> data(r,theta,phi)
        """
        self.M=int(self.N)
        self.volume_cellule = np.abs(self.nbc_to_pos(0)[0] - self.nbc_to_pos(1)[0]) **3 
        self.field_sph = np.zeros((self.M+1,self.M+1,self.M+1))
        self.table_nbc_to_rtp__nbc, self.table_nbc_to_rtp__rtp, self.table_nbc_to_irtp__irtp = [], [], []
        
        self.R_max = (self.pos_to_sph(self.nbc_to_pos(0))[0]) / mt.sqrt(3) # pos_to_sph(nb_to_pos(0)) : cellule dans un COIN. On veut sur la frontiere, mais atteingnable pour tout x,y,z 
        self.theta_max = np.pi
        self.phi_max = 2*np.pi
        
        self.dphi = self.phi_max/self.M
        self.dtheta = self.theta_max/self.M
        self.dr = self.R_max/(self.M)
        
        self.R = np.linspace(0,self.R_max,self.M+1)
        self.Theta = np.linspace(0,self.theta_max,self.M+1)
        self.Phi = np.linspace(0,self.phi_max,self.M+1)
        
        # for i in nb[int(N**3/2)-6:int(N**3/2)+6]:
        for i in self.nb:
            r, theta, phi = self.pos_to_sph(self.nbc_to_pos(i))
            if r <= self.R_max:
                diff_r, diff_theta, diff_phi = np.abs(self.R - r),np.abs(theta-self.Theta),np.abs(phi-self.Phi)
                
                indice_r = np.argwhere(diff_r == np.min(diff_r))[0][0]
                indice_theta = np.argwhere(diff_theta == np.min(diff_theta))[0][0]
                indice_phi = np.argwhere(diff_phi == np.min(diff_phi))[0][0]
                
                if rosette_board:
                   self.table_nbc_to_rtp__nbc += [int(i)] 
                   self.table_nbc_to_rtp__rtp += [(r,theta,phi)]
                   self.table_nbc_to_irtp__irtp += [(indice_r,indice_theta,indice_phi)]
                
                self.field_sph[indice_r,indice_theta,indice_phi] = self.data[int(i)]
