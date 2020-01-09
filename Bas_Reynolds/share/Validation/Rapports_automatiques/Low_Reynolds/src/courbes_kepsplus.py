#Script de recuperation des donnees u* et calcul de U+ et y+

import optparse
#import math
from math import *

import os

def properties():
    # ouverture des fichiers
    nomFic = 'propertiesGeometry.dat'

    fic = open(nomFic,'r')

    # lecture de ligne -> entetes
    fichier = fic.readlines()

    ligne = fichier [0]
    tLigne = ligne.split()
    mu=float(tLigne[0])
    rho=float(tLigne[1])
    h=float(tLigne[2])
    L=float(tLigne[3])

    fic.close()
    return mu,rho,h,L

#fonction de recuperation de la colonne u+
def utau():
    # ouverture des fichiers
    nomFicUstar = 'u_tau.dat'

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    u=float(tLigne[3])
    ficUstar.close()
    return u


def ecritureFichier(u,nomFicData,nomFicData2,rho,mu):
    #ecriture du fichier pour gnuplot
    nomFic = 'courbe_kepsplus.dat'
    ficRead = open(nomFicData,'r')
    ficRead2 = open(nomFicData2,'r')
    fichier = open(nomFic, 'w')

    #recuperation des valeurs. attention on suppose que ficRead et ficRead2 ont la meme structure
    fin=False
    while not fin:
        ligne = ficRead.readline()
        ligne2 = ficRead2.readline()
        if not ligne:
            fin=True
        tLigne = ligne.split()
        tLigne2 = ligne2.split()
        if len(tLigne)>0:
            y=float(tLigne[0])
            yPlus=y*rho*float(u)/mu
            kPlus=float(tLigne[1])/(u*u)
            epsPlus=mu*float(tLigne2[1])/(rho*u*u*u*u)
            Ret=kPlus*kPlus/(epsPlus+1e-15)
            fmu=exp(-3.4/((1+Ret/50)*(1+Ret/50)))
            fichier.write(' %18.8f  %18.8f %18.8f %18.8f\n' % (yPlus,kPlus,epsPlus,fmu))

    fichier.close()
    ficRead.close()
    ficRead2.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    mu,rho,h,L = properties()

    #recuperation des donnees
    u = utau()

    nomFic = args[0]
    nomFic2 = args[1]
    #ecriture du fichier gnuplot
    ecritureFichier(u,nomFic,nomFic2,rho,mu)
