# -*- coding: utf8

import os
import numpy
import subprocess
import scipy
from scipy import signal
import pandas
import glob

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import matplotlib.animation as animation


def pas_de_bulle_a_cheval(di):
    """ 
        di = dx(chi)
        Prend en entree un tableau [nsondes,ntemps,npoints]
        s'assure qu'on commence di par une valeur negative
        et qu'on termine di par une valeur positive
        ==> s'assure qu'aucune bulle ne traverse les frontieres.
            (liquide : 1; gaz : 0)
    """
    nsondes = di.shape[0]
    ntemps = di.shape[1]
    npoints = di.shape[2]
    
    inv_di = numpy.flip(di,axis=2)
    for s in range(nsondes):
        for t in range(ntemps):
            if numpy.argmax(di<-0.5)>numpy.argmax(di>0.5):
                di[s,t,numpy.argmax(di<-0.5)] = 0.
                
            if numpy.argmax(inv_di>0.5) > numpy.argmax(inv_di<-0.5):
                di[s,t,npoints-numpy.argmax(inv_di>0.5)-1] = 0.
    
    return di



def Gauss(abscisses,moyenne=0,ecart_type=1):
    arg = - (abscisses-moyenne)**2 / (ecart_type)**2
    pre_facteur = 1. / (numpy.sqrt(2*numpy.pi) * ecart_type)
    corps = pre_facteur*numpy.exp(arg)
    post_facteur = numpy.sum(corps)
    res = corps / post_facteur
    return(res)

def coordSonde(chemin_sonde):    
    sonde1 = chemin_sonde
    if os.path.isfile(sonde1):
        # On isole la ligne 2 avec les coordonnées
        subprocess.call(['sed -n "2{p;q}" '+sonde1+" > /tmp/fic.son"], shell=True)
        # Suppression du mot Temps
        subprocess.call(['sed -i "s/# Temps //g" /tmp/fic.son'], shell=True)
        # Passage à la ligne de chaque triplet de coordonnées et suppression de x=
        subprocess.call(['sed -i "s/x=/\\n/g" /tmp/fic.son'], shell=True)
        # Remplacement de y= par une virgule
        subprocess.call(['sed -i "s/ y= /,/g" /tmp/fic.son'], shell=True)
        # Remplacement des z= par une virgule
        subprocess.call(['sed -i "s/ z= /,/g" /tmp/fic.son'], shell=True)
        # Nettoyage des espaces restantes
        subprocess.call(['sed -i "s/ //g" /tmp/fic.son'], shell=True)
        # Ajout d'un titre par colonne pour pandas
        subprocess.call(['sed -i "1cx,y,z" /tmp/fic.son'], shell=True)
        
        return pandas.read_csv("/tmp/fic.son")
    else:
        print("Le fichier "+chemin_sonde+" n'existe pas.")
        
def coordSonde_without_init_and_end(chemin_sonde):    
    sonde1 = chemin_sonde
    if os.path.isfile(sonde1):
        # On isole la ligne 2 avec les coordonnées
        subprocess.call(['sed -n "2{p;q}" '+sonde1+" > /tmp/fic.son"], shell=True)
        # Suppression du mot Temps
        subprocess.call(['sed -i "s/# Temps //g" /tmp/fic.son'], shell=True)
        # Passage à la ligne de chaque triplet de coordonnées et suppression de x=
        subprocess.call(['sed -i "s/x=/\\n/g" /tmp/fic.son'], shell=True)
        # Remplacement de y= par une virgule
        subprocess.call(['sed -i "s/ y= /,/g" /tmp/fic.son'], shell=True)
        # Remplacement des z= par une virgule
        subprocess.call(['sed -i "s/ z= /,/g" /tmp/fic.son'], shell=True)
        # Nettoyage des espaces restantes
        subprocess.call(['sed -i "s/ //g" /tmp/fic.son'], shell=True)
        # Suppression de la premiere ligne
        subprocess.call(['sed -i \'1d;$d\' /tmp/fic.son'], shell=True)
        # Ajout d'un titre par colonne pour pandas
        subprocess.call(['sed -i "1cx,y,z" /tmp/fic.son'], shell=True)
        
        return pandas.read_csv("/tmp/fic.son")
    else:
        print("Le fichier "+chemin_sonde+" n'existe pas.")        
        
        

def loadSondeSeg(nom):
    print(" = Import de la sonde ",)
    print(nom+"...")
    return numpy.loadtxt(nom)

def loadPartSondeSeg(nom,it0,itN):
    print(" = Import de la sonde ",)
    print(nom+"...")
    return numpy.loadtxt(nom)[it0:itN,:]

def getNbPointNbTemps(sonde):
    # la première ligne est parfois pleine de 0 ou de 1
    if numpy.min(sonde[0, 1:]) == 0 and numpy.max(sonde[0, 1:]) == 0:
        sonde = numpy.delete(sonde, 0, 0)  
    # Le temps est à la première colonne
    temps = sonde[:, 0]
    # On retire la colonne du temps
    sonde = sonde[:, 1:]
    
    return(sonde.shape[1],sonde.shape[0])
   

class Donnees():
    def __init__(self):
        # Denombrement
        self.nsondes, self.npoints, self.ntemps = 2,2,2        
        # Tableaux de valeurs + dx
        self.brute = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.autocovariance = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre1 = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre2 = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre_welch = numpy.zeros((self.nsondes, self.ntemps, int(self.npoints/2+1)))
        self.deltaEspace = numpy.zeros((self.nsondes))
        # marqueur spécifiques (get_regular_fluctu_[time;space] et compute_correlation)
        self.add_T = True
        self.add_S = True
                        
    def initialize_datas(self):
        """ A appeler apres avoir appele get_the_n et broadcast_the_n de la classe Grandeur """
        self.brute = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.autocovariance = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre1 = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre2 = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre_welch = numpy.zeros((self.nsondes, self.ntemps, int(self.npoints/2+1)))
        self.frequences_welch = numpy.zeros((self.nsondes, self.ntemps, int(self.npoints/2+1)))
        self.deltaEspace = numpy.zeros((self.nsondes)) 
        
    def initialize_crude(self):
        """ A appeler apres avoir appele get_the_n et broadcast_the_n de la classe Grandeur """
        self.brute = numpy.zeros((self.nsondes, self.ntemps, self.npoints))        
        self.deltaEspace = numpy.zeros((self.nsondes)) 
         
        
    def initialize_data_after_get_crude(self):
        self.autocovariance = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre1 = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre2 = numpy.zeros((self.nsondes, self.ntemps, self.npoints))
        self.spectre_welch = numpy.zeros((self.nsondes, self.ntemps, int(self.npoints/2)+1))
        self.frequences_welch = numpy.zeros((self.nsondes, self.ntemps, int(self.npoints/2)+1))        
                      
    
    def get_whole_brute(self):
        """ Arrondi le relevé de la sonde """
        self.brute = numpy.rint(self.brute)
    
    def get_fancy_fluctu(self):
        """ definition alternative des fluctuations """
        self.fluctuations = (self.brute - numpy.mean(self.brute,axis=2)[...,numpy.newaxis]) / (1e-20+numpy.std(self.brute,axis=2)[...,numpy.newaxis])
        
    def get_all_stats(self):
        """ statistiques calculees sur le tableau totalemnt applatit en 1D """
        self.moyenne = self.brute.flatten().mean()
        self.variance = self.brute.flatten().var()
        self.skewness = scipy.stats.skew(self.brute.flatten())
        self.kurtosis = scipy.stats.kurtosis(self.brute.flatten())
        
    def get_all_stats_scipy(slef):
        """stats aura cette tete : (nobs=, minmax=( , ), mean=, variance=, skewness=, kurtosis="""
        self.stats = stats.describe(self.brute.flatten())
        
    def get_all_stats_time(self):
        """ statistiques calculees sur la dimension temporelle uniquement """
        self.moyenne_T = self.brute.reshape(self.nsondes*self.ntemps,self.npoints).mean(axis=0)
        self.variance_T = self.brute.reshape(self.nsondes*self.ntemps,self.npoints).var(axis=0)
        self.skewness_T = scipy.stats.skew(self.brute.reshape(self.nsondes*self.ntemps,self.npoints),axis=0)
        self.kurtosis_T = scipy.stats.kurtosis(self.brute.reshape(self.nsondes*self.ntemps,self.npoints),axis=0)
    
    def get_all_stats_space(self):
        """ statistiques calculees sur la dimension spatiale uniquement """
        self.moyenne_S = self.brute.transpose(1,0,2).reshape(self.ntemps,self.nsondes*self.npoints).mean(axis=1)
        self.variance_S = self.brute.transpose(1,0,2).reshape(self.ntemps,self.nsondes*self.npoints).var(axis=1)
        self.skewness_S = scipy.stats.skew(self.brute.transpose(1,0,2).reshape(self.ntemps,self.nsondes*self.npoints),axis=1)
        self.kurtosis_S = scipy.stats.kurtosis(self.brute.transpose(1,0,2).reshape(self.ntemps,self.nsondes*self.npoints),axis=1)
    
        
    def get_regular_fluctu(self):
        """ definition 'classique' des fluctuations """
        self.fluctuations = self.brute - self.moyenne
        
    def get_regular_fluctu_time(self):
        """ fluctuations par rapport à la moyenne temporelle, on ajoute un nouvel attribut à la méthode
            si l'attribut add_T est à True. Sinon, on écrit sur la variable fluctuations
        """
        if self.add_T : 
            self.fluctuations_T = self.brute - self.moyenne_T[numpy.newaxis,numpy.newaxis,:]
            self.fluctuations_T2 = self.brute - self.moyenne
        else : 
            self.fluctuations = self.brute - self.moyenne_T[numpy.newaxis,numpy.newaxis,:]
        
    def get_regular_fluctu_space(self):
        """ fluctuations par rapport à la moyenne spatiale, on ajoute un nouvel attribut à la méthode
            si l'attribut add_S est à True. Sinon, on écrit sur la variable fluctuations
        """
        if self.add_S :  
            self.fluctuations_S = self.brute - self.moyenne_S[numpy.newaxis,:,numpy.newaxis]
            self.fluctuations_S2 = self.brute - self.moyenne
        else : 
            self.fluctuations = self.brute - self.moyenne_S[numpy.newaxis,:,numpy.newaxis]
    

class Grandeur(Donnees):
    def __init__(self,Donnees):
        Donnees.__init__()
        # Informations generales
        repoSondes = "/volatile/RANDOM_TRIPERIO/a015/SONDES/CHI"
        formule = r'$\latex$'
        nom_grandeur = 'qu est ce que c est'
        nom_sonde = 'champ tel que dans le nom du .son'
        nom_fichier = 'le prefixe de tous les fichiers sonde'
        direction_sonde = 'S_X, S_Y ou S_Z'
        liste_nomSondes = []
        # Denombrement
        nsondes, npoints, ntemps = 2,2,2
        
    def get_nsondes(self):
        liste_sondes = glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")
        self.nsondes = len(liste_sondes)
    def get_ntemps(self):
        print(self.repoSondes+"/*"+self.direction_sonde+"*"+self.nom_sonde+".son")
        sonde = loadPartSondeSeg(glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")[-1],self.it_deb,self.it_fin)
        self.ntemps = sonde.shape[0]
    def get_npoints(self):
        sonde = loadPartSondeSeg(glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")[-1],self.it_deb,self.it_fin)
        self.npoints = sonde.shape[1]-1
        
    def get_the_n_old(self):
        """ Met a jour la valeur des n pour la Grandeur"""
        self.get_nsondes()
        self.get_ntemps()
        self.get_npoints()

    def get_the_n(self):
        """ Met a jour la valeur des n pour la Grandeur"""
        liste_sondes = glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")
        print('loadPartSondeSeg(glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")')
        print(glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son"))
        sonde = loadPartSondeSeg(glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")[-1],self.it_deb,self.it_fin)
        self.nsondes = len(liste_sondes)
        self.ntemps = int(sonde.shape[0])
        self.npoints = int(sonde.shape[1]-1)
    
                
    def get_crude_from_probe(self):
        """ Stocke la donnée de la sonde dans l'objet dedie et stocke les coordoonnees des points de releve """
        bruteT = numpy.zeros((self.nsondes,self.ntemps,self.npoints+1))
        liste_sondes = glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")
        for efic, fic in enumerate(liste_sondes[:]): 
            bruteT[efic,:,:] = loadPartSondeSeg(fic,self.it_deb,self.it_fin)
            nomSondes = fic.replace(self.repoSondes, "").replace("SONDES-X", "")
            nomSondes = nomSondes.replace("/", "").replace(self.nom_fichier+'_'+self.direction_sonde, "")
            nomSondes = nomSondes.replace('_'+self.nom_grandeur,"")
            nomSondes = nomSondes.replace('.son',"")
            self.liste_nomSondes.append(nomSondes)
                
            # ATTN il faudra le remettre
            # la première ligne est parfois pleine de 0 ou de 1
            # if numpy.min(bruteT[efic,0, 1:]) == 0 and numpy.max(bruteT[efic,0, 1:]) == 0:
                # bruteT = numpy.delete(bruteT[efic,:,:], 0, 0)  
            
            self.coord = coordSonde(fic)
            # Détermination de l'orientation de la sonde
            for i in self.coord.keys():
                if self.coord[i][2] - self.coord[i][1] != 0.:
                    self.deltaEspace[efic] = self.coord[i][10] - self.coord[i][9]
                    
        for efic, fic in enumerate(liste_sondes[:]):   
            # On retire la colonne du temps
            self.brute[efic,:,:] = bruteT[efic,:,1:]
    
    def get_crude_from_same_crude(self,same):
        """
        Dans le cas ou on a deja lu la sonde pour un autre objet,
        plutot que de re-parcourir le fichier,
        on copie la donnee directement
        """
        self.brute = same.brute.copy()
        self.coord = same.coord.copy()
        self.deltaEspace = same.deltaEspace
                    
    def get_rho_from_crude(self,rho_l,rho_v,indicatrice_l):
        self.formule = r'$\rho$'
        self.brute = indicatrice_l.brute*rho_l + (1-indicatrice_l.brute)*rho_v
        self.deltaEspace = indicatrice_l.deltaEspace
        self.npoints, self.ntemps, self.nsondes = indicatrice_l.npoints, indicatrice_l.ntemps, indicatrice_l.nsondes
        
    def get_rho_u_from_crude(self,rho,u):
        self.formule = r'$\rho u$'
        self.brute = rho.brute*u.brute
        self.deltaEspace = u.deltaEspace
        self.npoints, self.ntemps, self.nsondes = u.npoints, u.ntemps, u.nsondes
        
    def locate_crude_like(self,Grandeur):
        """ 
        rho est relevé aux centres des éléments.
        u_i est relevé sur les faces de normale i.
        Pour construire rho*u_i, il  faut relocaliser rho et u_i au meme point exactement
        /!\ : On ne peut relocaliser que dans la direction de la sonde !
              --> rho * u_y pour la sonde le long de x NE PEUT PAS etre reconstruit.
        """
        """
        Le premier relevé est fait n'importe comment. On le retire, ce qui force à perde non pas
        un mais 3 points.
        """
        """
        Interpolation linéaire : f(x) = (an-ao)/(xn-xo) * (x-xo) + ao
        """
        direction_Grandeur = Grandeur.direction_vitesse[-1].swapcase()
        print("direction_Grandeur : "+direction_Grandeur)
        print("Grandeur.direction_vitesse : "+Grandeur.direction_vitesse)
        print("self.direction_sonde : "+self.direction_sonde)
        # A changer : on prepare la derivee dans la direction de la sonde la
        coord_plus = numpy.array(self.coord[direction_Grandeur][:-2])
        coord_moins = numpy.array(self.coord[direction_Grandeur][1:-1])
        print(direction_Grandeur)
        # A changer : on prepare la derivee dans la direction de la sonde la
        dcoord = coord_plus - coord_moins
        print(dcoord)
        # A changer : on prepare la derivee dans la direction de la sonde la
        dbrute = self.brute[...,:-2] - self.brute[...,1:-1]
        dposition = (numpy.array(Grandeur.coord[direction_Grandeur][1:-1]) - coord_moins)[numpy.newaxis,numpy.newaxis,...]
        self.brute_reloc = dbrute/dcoord * (dposition) + self.brute[...,1:-1]
        print(":::::::: relocation :::::::::::::")
        print(":::::::: Tout est ok normalement :::::::::::::")
        
        # Cette partie rabote le premier et dernier point de la sonde car elles relevent un point absurde au début du segment
        # Il serai plus judicieux de l'accoller à la méthode "after_having_relocated"
        fic = glob.glob(Grandeur.repoSondes+"/*"+Grandeur.direction_sonde+"*"+"*"+Grandeur.identification+"*"+Grandeur.nom_sonde+"*"+".son")[0]
        self.coord = coordSonde_without_init_and_end(fic)
        self.npoints-=2
        
        
    def locate_crude_coaxial(self,Grandeur):
        debug = False
        """ 
        rho est relevé aux centres des éléments.
        u_i est relevé sur les faces de normale i.
        Pour construire rho*u_i, il  faut relocaliser rho et u_i au meme point exactement
        /!\ : On ne peut relocaliser que dans la direction de la sonde !
              --> rho * u_y pour la sonde le long de x NE PEUT PAS etre reconstruit.
        """
        """
        Le premier relevé est fait n'importe comment. On le retire, ce qui force à perde non pas
        un mais 3 points.
        """
        """
        Interpolation linéaire : f(x) = (an-ao)/(xn-xo) * (x-xo) + ao
        """
        direction_Grandeur = Grandeur.direction_vitesse[-1].swapcase()
        if debug : 
            print("direction_Grandeur : "+direction_Grandeur)
            print("Grandeur.direction_vitesse : "+Grandeur.direction_vitesse)
            print("self.direction_sonde : "+self.direction_sonde)
        # A changer : on prepare la derivee dans la direction de la sonde la
        coord_plus = numpy.array(self.coord[direction_Grandeur][:-2])
        coord_moins = numpy.array(self.coord[direction_Grandeur][1:-1])
        if debug : print(direction_Grandeur)
        # A changer : on prepare la derivee dans la direction de la sonde la
        dcoord = coord_plus - coord_moins
        if debug : print(dcoord)
        # A changer : on prepare la derivee dans la direction de la sonde la
        dbrute = self.brute[...,:-2] - self.brute[...,1:-1]
        dposition = (numpy.array(Grandeur.coord[direction_Grandeur][1:-1]) - coord_moins)[numpy.newaxis,numpy.newaxis,...]
        self.brute_reloc = dbrute/dcoord * (dposition) + self.brute[...,1:-1]
        print(":::::::: relocation :::::::::::::")
        print(":::::::: Tout est ok normalement :::::::::::::")
        
        """
        # Cette partie rabote le premier et dernier point de la sonde car elles relevent un point absurde au début du segment
        # Il serai plus judicieux de l'accoller à la méthode "after_having_relocated"
        fic = glob.glob(Grandeur.repoSondes+"/*"+Grandeur.direction_sonde+"*"+"*"+Grandeur.identification+"*"+Grandeur.nom_sonde+"*"+".son")[0]
        self.coord = coordSonde_without_init_and_end(fic)
        self.npoints-=2 
        """       
        
    def locate_crude_crossed(self, CenterLine, Grandeur):
        debug = False
        """
            On localise la grandeur self sur la Grandeur en entrée de paramètre. Pas l'inverse.
        """
        """ 
         - rho est relevé aux centres des éléments.
         - u_i est relevé sur les faces de normale i.
        Pour construire rho*u_i, il  faut relocaliser rho et u_i au meme point exactement
        Avec la donnée d'un segment sonde, on ne peut relocaliser que dans la direction de la sonde !
        Pour palier ce problème, genSondes_GR.py prends en compte un paramètre "mode"(="." ou "v" ou "+")
        qui place un faisceau de segmetns (en coin pou "v", en croix pour "+") et qui permet de relocaliser
        rho * u_y pour la sonde le long de x par exemple.
        """
        """
        Le premier relevé est fait n'importe comment. On le retire, ce qui force à perde non pas
        un mais 3 points.
        """
        """
        Interpolation linéaire : f(x) = (an-ao)/(xn-xo) * (x-xo) + ao
        """
        direction_Grandeur = Grandeur.direction_vitesse[-1].swapcase()
        if debug :
            print("direction_Grandeur : "+direction_Grandeur)
            print("Grandeur.direction_vitesse : "+Grandeur.direction_vitesse)
            print("self.direction_sonde : "+self.direction_sonde[2].swapcase())
        # A changer : on prepare la derivee dans la direction de la sonde la
        if direction_Grandeur==self.direction_sonde[2].swapcase() : 
            print("# Granderur.direction_vitesse et self.direction sont sont identiques : appelle plutot locate_crude_like")
            # print(direction_Grandeur, self.direction_sonde[2].swapcase())
        else :
            print("# Colocalisation croisée")
            # L'attribut "plus" ou "moins" n'est pas vraiment pertinent ici mais n'engendre pas de contre-sens pour autant
            # Seul l'ordre des différences (dcoord et dbrute) importe
            coord_plus = numpy.array(self.coord[direction_Grandeur][:])
            if debug:    
                print("coord_plus")
                print(self.direction_sonde,self.identification)
                print(coord_plus)
            coord_moins = numpy.array(CenterLine.coord[direction_Grandeur][:])
            if debug:    
                print("coord_moins")
                print(CenterLine.direction_sonde,CenterLine.identification)
                print(coord_moins)
            dcoord = coord_plus - coord_moins
            dbrute = self.brute[...,:] - CenterLine.brute[...,:]
            dposition = (numpy.array(Grandeur.coord[direction_Grandeur][:]) - coord_moins)[numpy.newaxis,numpy.newaxis,...]
            if debug : print(dcoord)
            self.brute = dbrute/dcoord * (dposition) + self.brute[...,:]
            print(":::::::: relocation :::::::::::::")
            print(":::::::: Tout est réglé normalement :::::::::::::")
        
        """
        Cette partie rabote le premier et dernier point de la sonde car elles relevent un point absurde au début du segment
        cette partie ne sert donc que pour les relocalisation "coaxiales" (vitesse dans la direction de la sonde)
        # Verifier si c'est self.identification ou Grandeur.identification qu'il faut ici en vrai..
        fic = glob.glob(Grandeur.repoSondes+"/*"+Grandeur.direction_sonde+"*"+"*"+Grandeur.identification+"*"+Grandeur.nom_sonde+"*"+".son")[0]
        self.coord = coordSonde_without_init_and_end(fic)
        self.npoints-=2
        """
        
        
    def locate_crude_all(self, CenterLine, Grandeur):
        direction_Grandeur = Grandeur.direction_vitesse[-1].swapcase()
        print("# direction_Grandeur : "+direction_Grandeur)
        print("# self.direction_sonde : "+self.direction_sonde)
        # A changer : on prepare la derivee dans la direction de la sonde la
        if direction_Grandeur==self.direction_sonde[2].swapcase() : 
            print("# Locate coax ")
            self.locate_crude_coaxial(Grandeur)
        else :
            print("# Locate crossed")
            self.locate_crude_crossed(CenterLine, Grandeur)
        self.coord_reloc = Grandeur.coord
            
        
        
    def after_having_relocated(self):
        """
        Pour raboter la grandeur qui a servi de modele pour relocaliser
        """
        direction = self.direction_sonde[2].swapcase()
        self.coord[direction] = self.coord[direction][1:-1]
        self.npoints-=2
        self.brute = self.brute[...,1:-1]
        fic = glob.glob(self.repoSondes+"/*"+self.direction_sonde+"*"+"*"+self.identification+"*"+self.nom_sonde+".son")[0]
        self.coord = coordSonde_without_init_and_end(fic)        
        
        
def compute_correlations(donnee):
    """
    Calcule uniquement les correlations de la donnee.
          - Correlations temporelles, avec les moyennes et variances temporelles
          - Correlations spatiales, avec les moyennes et variances spatiales uniquement
          => Ce choix est complètement discutable sur la pertinance physique.
             Il est fait surtout parce qu'il permet de croiser l'axe des abscisses à un moment
             (ce qui n'est pas nécessairement le cas si on choisi la moyenne de TOUTE la série statistique)
    """        
    if not(donnee.add_S) : 
        donnee.fluctuations_S = donnee.fluctuations
        donnee.fluctuations_S2 = donnee.fluctuations
        print("Attention :  est-ce que la fluctuation a ete construite avec get_regular_fluctu_space ?\n\
        si non, la correlation spatiale a toutes les chances de ne rien vouloir dire...")

    donnee.correlation_spatiale = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    donnee.correlation_spatiale2 = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    for efic in range(donnee.nsondes):
        for t in range(1,donnee.ntemps):
            # FT(<fluctu fluctu>)
            donnee.correlation_spatiale[efic,t,:] = numpy.correlate(donnee.fluctuations_S[efic,t,:], donnee.fluctuations_S[efic,t,:], mode='same')
            donnee.correlation_spatiale2[efic,t,:] = numpy.correlate(donnee.fluctuations_S2[efic,t,:], donnee.fluctuations_S2[efic,t,:], mode='same')
            
        
    if not(donnee.add_T) : 
        donnee.fluctuations_T =  donnee.fluctuations
        donnee.fluctuations_T2 = donnee.fluctuations
        print("Attention :  est-ce que la fluctuation a ete construite avec get_regular_fluctu_time ?\n\
        si NON, la correlation temporelle a toutes les chances de ne rien vouloir dire...")    

    donnee.correlation_temporelle = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    donnee.correlation_temporelle2 = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    for efic in range(donnee.nsondes):
        for x in range(1,donnee.npoints):
            # Dans le cas de donnees tres bruitees on peut vouloir moyenner notre fluctuation
            #    Avec des produits de ocnvolutions bien choisis (scipy.fftconvolve) on pourrait directement faire l'opération sur toutes les dimentions je pense
            donnee.fluctuations_T[efic,:,x] =  numpy.convolve(donnee.fluctuations_T[efic,:,x],donnee.smoothing_window/(float(len(donnee.smoothing_window))),mode='same')
            donnee.fluctuations_T2[efic,:,x] = numpy.convolve(donnee.fluctuations_T2[efic,:,x],donnee.smoothing_window/(float(len(donnee.smoothing_window))),mode='same')
            
            donnee.correlation_temporelle[efic,:,x] = numpy.correlate(donnee.fluctuations_T[efic,:,x], donnee.fluctuations_T[efic,:,x], mode='same')
            donnee.correlation_temporelle2[efic,:,x] = numpy.correlate(donnee.fluctuations_T2[efic,:,x], donnee.fluctuations_T2[efic,:,x], mode='same')
    
    # Moyennes
    donnee.correlation_spatiale = numpy.mean(numpy.mean(donnee.correlation_spatiale,axis=1),axis=0)
    donnee.correlation_spatiale2 = numpy.mean(numpy.mean(donnee.correlation_spatiale2,axis=1),axis=0)
    donnee.correlation_temporelle = numpy.mean(numpy.mean(donnee.correlation_temporelle,axis=2),axis=0)
    donnee.correlation_temporelle2 = numpy.mean(numpy.mean(donnee.correlation_temporelle2,axis=2),axis=0)
    
        
def compute_3spectra(donnee):
    """ ATTENTION : l'argument a entrer est une 'Grandeur'...
        Calcule les spectres 
            1 : TF(<fluctu fluctu>)
            2 : TF(fluctu).TF*(fluctu)
            Welch : 
    """
    fft_fluctuation = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    for efic in range(donnee.nsondes):
        for t in range(1,donnee.ntemps):
            # FT(<fluctu fluctu>)
            donnee.autocovariance[efic,t,:] = numpy.correlate(donnee.fluctuations[efic,t,:], donnee.fluctuations[efic,t,:], mode='same')
            donnee.spectre1[efic,t,:] = scipy.fft(donnee.autocovariance[efic,t,:])
            donnee.spectre1[efic,t,:] = numpy.absolute(donnee.spectre1[efic,t,:]) / (float((donnee.npoints**2)))
            # FT(fluctu) FT*(fluctu)
            fft_fluctuation[efic,t,:] = scipy.fft(donnee.fluctuations[efic,t,:])
            donnee.spectre2[efic,t,:] = fft_fluctuation[efic,t,:]*numpy.conjugate(fft_fluctuation[efic,t,:])
            donnee.spectre2[efic,t,:] = numpy.absolute(donnee.spectre2[efic,t,:]) / float((donnee.npoints**2))
            # Welch
            # donnee.spectre_welch[efic,t,:] = scipy.signal.welch(donnee.fluctuations[efic,t,:], 1./(donnee.deltaEspace[efic]),scaling='density')[1]
            # donnee.frequences_welch[efic,t,:] = scipy.signal.welch(donnee.fluctuations[efic,t,:], 1./(donnee.deltaEspace[efic]),scaling='density')[0]
    # Moyennes
    donnee.spectre1_moyen = numpy.mean(numpy.mean(donnee.spectre1,axis=1),axis=0)/donnee.variance
    donnee.spectre2_moyen = numpy.mean(numpy.mean(donnee.spectre2,axis=1),axis=0)/donnee.variance
    donnee.frequences_spectre = numpy.fft.fftfreq(donnee.npoints,numpy.mean(donnee.deltaEspace)).real
    # donnee.spectre_welch_moyen = numpy.mean(numpy.mean(numpy.mean(donnee.spectre_welch,axis=2),axis=1),axis=0)/donnee.variance
    # donnee.autocovariance_moyen = numpy.mean(numpy.mean(numpy.mean(donnee.autocovariance,axis=2),axis=1),axis=0)

def compute_spectrum1_AB(donneeA,donneeB):
    """ ATTENTION : l'argument a entrer est une 'Grandeur'...
    Calcule le spectre de A.B de la facon suivante :
    Sp(A.B) = TF[Corr(A;B)]
    """
    spectre1_moyen = numpy.zeros((donneeA.nsondes,donneeA.ntemps,donneeA.npoints))
    for efic in range(donneeA.nsondes):
        for t in range(1,donneeA.ntemps):
            covariance = numpy.correlate(donneeA.fluctuations[efic,t,:], donneeB.fluctuations[efic,t,:], mode='same')
            spectre1_moyen[efic,t,:] = scipy.fft(covariance)
            spectre1_moyen[efic,t,:] = numpy.absolute(spectre1_moyen[efic,t,:]) / (float((donneeA.npoints**2)))
                        
    variance =  (donneeA.brute*donneeB.brute).flatten().var()
    spectre1_moyen = numpy.mean(numpy.mean(spectre1_moyen,axis=1),axis=0) / variance
    frequences_spectre = numpy.fft.fftfreq(donneeA.npoints,numpy.mean(donneeA.deltaEspace)).real
    return(spectre1_moyen,frequences_spectre)
            
def compute_spectrum2_AB(donneeA,donneeB):
    """ ATTENTION : l'argument a entrer est une 'Grandeur'...
    Calcule le spectre de A.B de la facon suivante :
    Sp(A.B) = TF[A] . TF*[B]
      --> Il faudrai ajouter TF*[A] . TF[B] je pense...
    """
    spectre2_moyen = numpy.zeros((donneeA.nsondes,donneeA.ntemps,donneeA.npoints))
    for efic in range(donneeA.nsondes):
        for t in range(1,donneeA.ntemps):
            fftA = scipy.fft(donneeA.fluctuations[efic,t,:])
            fftB = scipy.fft(donneeB.fluctuations[efic,t,:])
            somme = fftA*numpy.conjugate(fftB) + fftB*numpy.conjugate(fftA)
            spectre2_moyen[efic,t,:] = somme
            # spectre2_moyen[efic,t,:] +=
            spectre2_moyen[efic,t,:] = numpy.absolute(spectre2_moyen[efic,t,:]) / (float((donneeA.npoints**2)))
                        
    variance =  (donneeA.brute*donneeB.brute).flatten().var()
    spectre2_moyen = numpy.mean(numpy.mean(spectre2_moyen,axis=1),axis=0) / variance
    frequences_spectre = numpy.fft.fftfreq(donneeA.npoints,numpy.mean(donneeA.deltaEspace)).real
    return(spectre2_moyen,frequences_spectre)
    
def compute_fourier_Aj_djBi(donneeA,donneeB):
    """ ATTENTION : l'argument a entrer est une 'Grandeur'...
    Calcule la grandeur TF(A.dB) de la facon suivante :
    TF(A.B) = TF(A) [X] TF(dB) = TF(A) [X] {k.TF(B)}
      --> Pratique quand A ou B comporte une derivee
    """
    # fftA = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    # fftB = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    fft_AdB = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    frequences_spectre = numpy.fft.fftfreq(donneeA.npoints,numpy.mean(donneeA.deltaEspace)).real
    for efic in range(donnee.nsondes):
        for t in range(1,donnee.ntemps):
            fft_A = scipy.fft(donneeA.fluctuations[efic,t,:])
            fftdB = scipy.fft(donneeB.fluctuations[efic,t,:])
            fftdB = frequences_spectre*fftdB
            fft_AdB[efic,t,:]  =                 scipy.signal.fftconvolve(fft_A,fftdB, mode='same')
            fft_AdB[efic,t,:] *= numpy.conjugate(scipy.signal.fftconvolve(fft_A,fftdB), mode='same')
            fft_AdB[efic,t,:] = numpy.absolute(fft_AdB[efic,t,:])
            fft_AdB[efic,t,:] /= (float((donneeA.npoints**2)))
            
    return(fft_AdB)
    
def compute_F_Ai__Fc_BjdjCi(donneeA,donneeB,donneeC):
    """ ATTENTION : l'argument a entrer est une 'Grandeur'...
    Calcule la grandeur TF(A.dB) de la facon suivante :
    TF(A.B) = TF(A) [X] TF(dB) = TF(A) [X] {k.TF(B)}
      --> Pratique quand A ou B comporte une derivee
      --> On recupere k selon la direction e la sonde de B.
          * C'est juste dans notre cas simple.
          * C'est faux si nos mailles ou points de releve des sondes
          ne sont pas uniformes.
    """
    # fftA = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    # fftB = numpy.zeros((donnee.nsondes,donnee.ntemps,donnee.npoints))
    fft_ABdC = numpy.zeros((donneeA.nsondes,donneeA.ntemps,donneeA.npoints))
    frequences_spectre = numpy.fft.fftfreq(donneeB.npoints,numpy.mean(donneeB.deltaEspace)).real
    for efic in range(donneeA.nsondes):
        for t in range(1,donneeA.ntemps):
            fft_A = scipy.fft(donneeA.fluctuations[efic,t,:])
            fft_B = scipy.fft(donneeB.fluctuations[efic,t,:])
            fft_dC = scipy.fft(donneeC.fluctuations[efic,t,:])
            fft_dC = frequences_spectre*fft_dC
            fft_BdC = scipy.signal.fftconvolve(fft_B,fft_dC, mode='same')
            
            fft_ABdC[efic,t,:] = fft_A*numpy.conjugate(fft_BdC)
            fft_ABdC[efic,t,:] = numpy.absolute(fft_ABdC[efic,t,:])
    # variance =  (donneeA.brute * donneeB.brute).flatten().var()
    # spectre3_moyen = numpy.mean(numpy.mean(spectre3_moyen,axis=1),axis=0)
    return(fft_ABdC)  
    
    
def compute_F_Ai__Fc_Bi(donneeA,donneeB): 
    """ ATTENTION : l'argument a entrer est une 'Grandeur'...
    calcule la grandeur TF(A).TF*(B)
    """
    fft_AB = numpy.zeros((donneeA.nsondes,donneeA.ntemps,donneeA.npoints))
    for efic in range(donneeA.nsondes):
        for t in range(1,donneeA.ntemps):
            fft_A = scipy.fft(donneeA.fluctuations[efic,t,:])
            fft_B = scipy.fft(donneeB.fluctuations[efic,t,:])
            
            fft_AB[efic,t,:] = fft_A*numpy.conjugate(fft_B)
            fft_AB[efic,t,:] = numpy.absolute(fft_AB[efic,t,:])
            
    return(fft_AB)
                                

class Vecteur():
    def __init__(self):
        # Informations generales
        direction_sonde = 'S_X, S_Y ou S_Z'
        dimensions = 'si on a u_x,u_y,u_z ou u_x,u_y ou ...'
        GrandeurX = 1
        GrandeurY = 1
        GrandeurZ = 1
        
    def scalar(self,Vecteur2):
        if self.dimensions == Vecteur2.dimensions:
            res = self.GrandeurX*Vecteur2.GrandeurX
        
            
            
            
