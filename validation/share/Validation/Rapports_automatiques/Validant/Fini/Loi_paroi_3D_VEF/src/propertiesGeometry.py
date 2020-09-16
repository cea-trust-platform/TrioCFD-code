#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties(nomFic):
    # ouverture des fichiers
    fic = open(nomFic,'r')

    chaines = ["Longueurs",
               "mu Champ_Uniforme",
               "rho Champ_Uniforme",] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="Longueurs":
                    Lx=float(tLigne[-3])
                    Ly=float(tLigne[-2])
                    Lz=float(tLigne[-1])
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[-1])
                if chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[-1])
    U=0.029
    fic.close()
    return Lx,Ly,Lz,mu,rho,U


def ecritureFichier(Lx,Ly,Lz,mu,rho,U):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.8f\n' % (Lx,Ly,Lz,mu,rho,U))
    fichier.close()

if __name__ == '__main__':

    #recuperation du fichier data
    #derniere ligne du ls
    ficLS = os.popen('ls *.data')
    lignes = ficLS.readlines()
    Ligne = lignes[0]
    #suppression du \n en fin de nom
    nomFic = Ligne[:len(Ligne)-1]
    Lx,Ly,Lz,mu,rho,U = properties(nomFic)

    #ecriture du fichier
    ecritureFichier(Lx,Ly,Lz,mu,rho,U)
