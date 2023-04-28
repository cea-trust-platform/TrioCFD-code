#Script de recuperation des proprietes physiques et de la geometrie

import os
import optparse
import sys
import math

def readPropertiesData():
    # ouverture des fichiers
    nomFic = 'test.data'

    fic = open(nomFic,'r')

    chaines = ["Longueurs",
               "mu Champ_Uniforme",
               "rho Champ_Uniforme",
               "# Vitesse_entree_calcul"] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="Longueurs":
                    h=float(tLigne[2])
                    L=float(tLigne[1])
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[3])
                if chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[3])
                if chaine=="# Vitesse_entree_calcul":
                    U0=float(tLigne[2])
    fic.close()
    return mu,rho,h,L,U0
    return properties


def ecritureFichier(mu,rho,h,L,U0):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f\n' % (mu,rho,h,L,U0))
    fichier.close()

if __name__ == '__main__':


    mu,rho,h,L,U0= readPropertiesData()

    #ecriture du fichier
    ecritureFichier(mu,rho,h,L,U0)
