#Script de recuperation des donnees u* et calcul de U+ et y+

import optparse
import math
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

#fonction de recuperation de la colonne u*
def uStar():
    # ouverture des fichiers
    nomFicUstar ='Prepare_pb_Ustar.face'

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    u=float(tLigne[5])
    ficUstar.close()

    return u

def ecrire_profils_vitesse(sonde,u_tau,mu,rho):

    ficRead = open(sonde, 'r')

    nomFic = sonde+'.dat'
    ficWrite = open(nomFic, 'w')

    nomFic2 = sonde+'_dimensionless.dat'
    ficWrite2 = open(nomFic2, 'w')

    ligne = ficRead.readlines()
    n=len(ligne)

    for i in range (1,n):

        tligne=ligne[i].split()
        ficWrite.write('%18.8f %18.8f\n' % (float(tligne[0]),float(tligne[1])))
        ficWrite2.write('%18.8f %18.8f\n' % (float(tligne[0])*u_tau*rho/mu,float(tligne[1])/u_tau))

    ficWrite.close()
    ficWrite2.close()
    ficRead.close()


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    mu,rho,h,L = properties()

    #recuperation des donnees
    u = uStar()

    #ecriture du profils de vitesse Vx=f(y) et U+=f(y+)
    fichierSonde='Prepare_SONDE_VIT.coupe'
    ecrire_profils_vitesse(fichierSonde,u,mu,rho)
