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
def uStar(dim):
    # ouverture des fichiers
    nomFicUstar = dim+'_keps_pb_Ustar.face'

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


def ecritureFichier(dim,u,nomFicData,rho,mu):
    #ecriture du fichier pour gnuplot
    nomFic = 'courbe_reichardt.dat'
    if dim=="3D":
        Re=100000

    ficRead = open(nomFicData,'r')
    fichier = open(nomFic, 'w')

    # lecture de ligne -> entetes
    i=0
    while i<17:
        ligne = ficRead.readline()
        i=i+1

    #recuperation des valeurs
    fin=False
    while not fin:
        ligne = ficRead.readline()
        if not ligne:
            fin=True
        tLigne = ligne.split()
        if len(tLigne)>0:
            y=float(tLigne[0])
            yPlus=y*rho*float(u)/mu
            uPlus=float(tLigne[1])/u
            fichier.write(' %18.8f %18.8f\n' % (yPlus,uPlus))

    fichier.close()
    ficRead.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    mu,rho,h,L = properties()

    #recuperation des donnees
    u = uStar(args[0])

    #recuperation du fichier avec le dernier temps dans le nom
    #derniere ligne du ls
    ficLS = os.popen('ls -rt '+args[1]+'_*')
    lignes = ficLS.readlines()
    derLigne = lignes[-1]
    #suppression du \n en fin de nom
    nomFic = derLigne[:len(derLigne)-1]

    #ecriture du fichier gnuplot
    ecritureFichier(args[0],u,nomFic,rho,mu)
