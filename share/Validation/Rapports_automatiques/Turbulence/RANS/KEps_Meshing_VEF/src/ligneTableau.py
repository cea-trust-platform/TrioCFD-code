#Script de recuperation des donnees y+, Re_tau et Relative error
import os

import optparse
import math

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

#fonction de recuperation de la colonne y+
def yPlus(dim):
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
    y=float(tLigne[3])
    ficUstar.close()

    return y

def retau(dim,rho,h,mu):
    # ouverture des fichiers
    nomFicRetau = 'reynolds_tau.dat'

    ficRetau = open(nomFicRetau,'r')

    # lecture de ligne -> entetes
    fichier = ficRetau.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    ReTau=float(tLigne[2])
    ficRetau.close()

    Re=rho*10*h/mu/2
    Puissance = math.pow(float(Re),(7.0/8.0))
    ReTheo=0.1754*Puissance

    return ReTau,ReTheo

def perte_charge(mu ):
    # ouverture des fichiers
    nomFic = os.popen('ls -rt *Pressure_Gradient_pb_*').readlines()[-1]
    nomFic = nomFic[:len(nomFic)-1] # Suppress /n

    fic = open(nomFic,'r')

    # lecture de ligne -> entetes
    fichier = fic.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    Perte=float(tLigne[1])
    fic.close()
    Re=100000
    Puissance = math.pow(float(Re),(7.0/8.0))
    Pertetheo=(0.1754*Puissance*mu)*(0.1754*Puissance*mu)
    return Perte,Pertetheo





def ecritureFichier(dim,y,Retau,RelError,Retheo,Perte,RelErrorPerte):
    #ecriture du fichier pour gnuplot
    nomFic = 'ligneTableau.dat'
    fichier = open(nomFic, 'w')
    fichier.write('_keps'+' %18.1f %18.1f %18.1f %18.1f %18.3f %18.1f\n' % (y,Retau,RelError,Retheo,Perte,RelErrorPerte))
    fichier.close()


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    mu,rho,h,L = properties()

    #recuperation des donnees
    y = yPlus(args[0])
    Retau,Retheo=retau(args[0],rho,h,mu)
    RelError=(abs(Retau-Retheo)/Retheo)*100
    Perte,Pertetheo=perte_charge(mu)
    RelErrorPerte=(abs(Pertetheo-Perte)/Pertetheo)*100
    #ecriture du fichier gnuplot
    ecritureFichier(args[0],y,Retau,RelError,Retheo,Perte,RelErrorPerte)
