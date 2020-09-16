#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties(nomFic):
    # ouverture des fichiers
    fic = open(nomFic,'r')

    chaines = ["Longueurs",
               "mu Champ_Uniforme",
               "rho Champ_Uniforme","extrudebord"] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="extrudebord":
                    h=2.
                    L=0.8
                if chaine=="Longueurs":
                    h=float(tLigne[2])
                    L=float(tLigne[1])
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[3])
                if chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[3])
    fic.close()
    return mu,rho,h,L


def ecritureFichier(mu,rho,h,L):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f\n' % (mu,rho,h,L))
    fichier.close()

if __name__ == '__main__':

    #recuperation du fichier data
    #derniere ligne du ls
    ficLS = os.popen('ls *.data')
    lignes = ficLS.readlines()
    derLigne = lignes[-1]
    #suppression du \n en fin de nom
    nomFic = derLigne[:len(derLigne)-1]
    mu,rho,h,L = properties(nomFic)

    #ecriture du fichier
    ecritureFichier(mu,rho,h,L)
