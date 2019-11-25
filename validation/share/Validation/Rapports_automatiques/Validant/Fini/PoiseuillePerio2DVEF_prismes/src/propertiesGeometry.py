#Script de recuperation des propretes physiques et de la geometrie et de la vitesse de base

import optparse
import math

def properties():
    # ouverture des fichiers
    nomFic = '../test.data'

    fic = open(nomFic,'r')

    chaines = ["# Longueurs",
               "mu Champ_Uniforme",
               "rho Champ_Uniforme",
               "vitesse Champ_fonc_xyz"] # Texte a rechercher
    print '*2'
    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="# Longueurs":
                    h1=float(tLigne[3])
                    print 'h1:%f' %h1
                    L1=float(tLigne[2])
                    print 'L1:%f' %L1
                if chaine=="mu Champ_Uniforme":
                    mu1=float(tLigne[3])
                    print 'mu1:%f' %mu1
                if chaine=="rho Champ_Uniforme":
                    rho1=float(tLigne[3])
                    print 'rho1:%f' %rho1
                if chaine=="vitesse Champ_fonc_xyz":
                    U01=float(tLigne[4])
                    print 'U01:%f' %U01
    print 'indice'
    fic.close()
    return mu1,rho1,h1,L1,U01


def ecritureFichier(mu,rho,h,L,U0):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f\n' % (mu,rho,h,L,U0))
    fichier.close()

if __name__ == '__main__':


    mu,rho,h,L,U0 = properties()

    print 'mu:%f' %mu

    print '*1'
    #ecriture du fichier
    ecritureFichier(mu,rho,h,L,U0)
