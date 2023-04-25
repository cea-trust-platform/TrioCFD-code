#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties(nomFic):
    # ouverture des fichiers
    fic = open(nomFic,'r')

    chaines = ["paroi_temperature_imposee Champ_front_uniforme",
               "T_ext Champ_Front_Uniforme",
               "mu Champ_Uniforme",
               "lambda Champ_Uniforme",
               "Longueurs",
               "Cp Champ_Uniforme",
               "rho Champ_Uniforme",
               "diam_hydr",
               "vitesse Champ_Uniforme"] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="Longueurs":
                    H=float(tLigne[-2])
                    L=float(tLigne[-1])
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[-1])
                if chaine=="lambda Champ_Uniforme":
                    lamb=float(tLigne[-1])
                if chaine=="Cp Champ_Uniforme":
                    Cp=float(tLigne[-1])
                if chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[-1])
                if chaine=="paroi_temperature_imposee Champ_front_uniforme" or chaine=="T_ext Champ_Front_Uniforme":
                    Tw=float(tLigne[-1])
                if chaine=="diam_hydr":
                    D=float(tLigne[-3])
                if chaine=="vitesse Champ_Uniforme":
                    u=abs(float(tLigne[-1]))
    fic.close()
    return H,L,mu,lamb,Cp,rho,Tw,D,u


def ecritureFichier(H,L,mu,lamb,Cp,rho,Tw,D,u):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f\n' % (H,L,mu,lamb,Cp,rho,Tw,D,u))
    fichier.close()

if __name__ == '__main__':

    #recuperation du fichier data
    #derniere ligne du ls
    ficLS = os.popen('ls *.data')
    lignes = ficLS.readlines()
    Ligne = lignes[0]
    #suppression du \n en fin de nom
    nomFic = Ligne[:len(Ligne)-1]
    H,L,mu,lamb,Cp,rho,Tw,D,u = properties(nomFic)

    #ecriture du fichier
    ecritureFichier(H,L,mu,lamb,Cp,rho,Tw,D,u)
