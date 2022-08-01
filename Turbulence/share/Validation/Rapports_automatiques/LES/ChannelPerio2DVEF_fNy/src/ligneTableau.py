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
    U0=float(tLigne[4])

    fic.close()
    return mu,rho,h,L,U0

#fonction de recuperation de la colonne y+
def yPlus():
    # ouverture des fichiers

    nomFicUstar = '2D_pb_Ustar.face'


    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    y=float(tLigne[3])
    utau_trio=float(tLigne[5])
    ficUstar.close()

    return y,utau_trio

#CEtte fonction calcul  aussi u_tau theorique a partir de ReTheo
def retau(rho,h,mu):
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
    ReTau=float(tLigne[3])
    ficRetau.close()

    Re=rho*10*h/mu/2
    Puissance = math.pow(float(Re),(7.0/8.0))
    ReTheo=0.1754*Puissance
    utau_theo = mu*ReTheo/rho/(h/2)

    return ReTau,ReTheo,utau_theo


def ecritureFichier(Ny,Reb,y,Retau,RelError,Retheo,utau_trio,utau_theo,u_tau_error,Pth,Ptr,Perror):
    #ecriture du fichier pour gnuplot
    nomFic = 'ligneTableau.dat'

    fichier = open(nomFic, 'w')
    fichier.write('%18.4f %18.4f %18.4f %18.4f %18.4f %18.4f %18.4f %18.4f %18.2f %18.4f %18.4f %18.2f\n' % (Ny,Reb,y,Retau,RelError,Retheo,utau_trio,utau_theo,u_tau_error,Pth,Ptr,Perror))
    fichier.close()


def Ptheo(utau_theo,L,h):
    Pth = float(utau_theo*utau_theo*L/(h/2))
    return Pth

def Ptrio(L):
    # ouverture des fichiers
    nomFic = os.popen('ls -rt *Pressure_Gradient_pb_periox').readlines()[-1]
    nomFic = nomFic[:len(nomFic)-1] # Suppress /n

    fic = open(nomFic,'r')

    # lecture de ligne -> entetes
    fichier = fic.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    ptrio=float(tLigne[1])*float(L)
    fic.close()

    return ptrio



if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    mu,rho,h,L,U0 = properties()

    #recuperation des donnees
    Ny=float(args[0])
    Reb=float(args[1])
    y,utau_trio = yPlus()
    Retau,Retheo,utau_theo=retau(rho,h,mu)
    RelError=(abs(Retau-Retheo)/Retheo)*100
    utau_error=float(abs((utau_theo-utau_trio)/utau_theo))*100

    Pth=Ptheo(utau_theo,L,h)
    Ptr=Ptrio(L)
    Perror=float(abs((Pth-Ptr)/Pth))*100



    #ecriture du fichier gnuplot
    ecritureFichier(Ny,Reb,y,Retau,RelError,Retheo,utau_trio,utau_theo,utau_error,Pth,Ptr,Perror)
