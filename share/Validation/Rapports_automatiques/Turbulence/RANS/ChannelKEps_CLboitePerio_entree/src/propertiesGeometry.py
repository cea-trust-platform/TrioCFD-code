#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties(nomFic):
    # ouverture des fichiers
    fic = open(nomFic,'r')

    chaines = ["Longueurs",
               "mu Champ_Uniforme",
               "rho Champ_Uniforme",
               "# Tailles",
               "vitesse Champ_fonc_xyz"] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="Longueurs":
                    h=float(tLigne[2])
                    L=float(tLigne[1])
                if chaine=="# Tailles":
                    h=float(tLigne[3])
                    L=float(tLigne[2])
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[3])
                if chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[3])
                if chaine=="vitesse Champ_fonc_xyz":
                    Ub=float(tLigne[4])
    fic.close()
    return mu,rho,h,L,Ub


def ecritureFichier(mu,rho,h,L,Ub):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f\n' % (mu,rho,h,L,Ub))
    fichier.close()

if __name__ == '__main__':

    #recuperation du fichier data
    #derniere ligne du ls
    ficLS = os.popen('ls *.data')
    lignes = ficLS.readlines()
    derLigne = lignes[-1]
    #suppression du \n en fin de nom
    nomFic = derLigne[:len(derLigne)-1]
    mu,rho,h,L,Ub = properties(nomFic)

    #ecriture du fichier
    ecritureFichier(mu,rho,h,L,Ub)
