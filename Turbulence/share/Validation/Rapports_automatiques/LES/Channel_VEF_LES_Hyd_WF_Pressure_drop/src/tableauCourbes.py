#Script de construction des tableaux Euler_Explicit et Crank_Nicholson et leurs courbes
import os

import optparse
import math

#fonction de recuperation de la colonne y+
def yPlus(scheme,U):
    # ouverture des fichiers
    nomFicUstar = 'canal_turbu_'+scheme+'_pb1_Ustar.face'

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

def Re(mu,rho,h,U):
    reh = float(rho)*float(U)*2*float(h)/float(mu)
    return reh

def Utheo(reh,mu,rho,h):
    rebulk = float(reh/4)
    retau = 0.1754*math.pow(float(rebulk),(7.0/8.0))
    utau = mu*retau/rho/(h/2)
    return utau

def Utrio(scheme,U):
    # ouverture des fichiers
    nomFicUstar = 'canal_turbu_'+scheme+'_pb1_Ustar.face'

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    utau=float(tLigne[5])
    ficUstar.close()

    return utau

def Uerr(Utheo,Utrio):
    Uerror = float(abs(Utheo-Utrio)/Utheo)*100
    return Uerror

def Ptheo(Utheo,L,h):
    Pth = float(Utheo*Utheo*L/(h/2))
    return Pth

def Ptrio(scheme,U,L):
    # ouverture des fichiers
    nomFic = os.popen('ls -rt *Pressure_Gradient_pb1_Entree').readlines()[-1]
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

def Perr(Ptheo,Ptrio):
    Uerror = float(abs(Ptheo-Ptrio)/Ptheo)*100
    return Uerror

def ecritureFichier(scheme,U,Re,y,Utheo,Utrio,UError,Ptheo,Ptrio,PError):
    #ecriture du fichier pour la ligne du tableau
    nomFicLigne = 'ligneTableau.dat'
    fichier = open(nomFicLigne, 'w')
    fichier.write('%s %18.2f %18.2f %e %e %18.2f %e %e %18.2f\n' % (U,Re,y,Utheo,Utrio,UError,Ptheo,Ptrio,PError))
    fichier.close()

    #ecriture du fichier pour la courbe dPtheo
    nomFicCourbedPtheo = '../CourbedPtheo_'+scheme+'.dat'
    fichier = open(nomFicCourbedPtheo, 'a')
    fichier.write('%18.8f %18.8f\n' % (Re,Ptheo))
    fichier.close()

    #ecriture du fichier pour la courbe dPtrio
    nomFicCourbedPtrio = '../CourbedPtrio_'+scheme+'.dat'
    fichier = open(nomFicCourbedPtrio, 'a')
    fichier.write('%18.8f %18.8f\n' % (Re,Ptrio))
    fichier.close()

    #ecriture du fichier pour la courbe Utau theo
    nomFicCourbeUtautheo = '../CourbeUtautheo_'+scheme+'.dat'
    fichier = open(nomFicCourbeUtautheo, 'a')
    fichier.write('%18.8f %18.8f\n' % (Re,Utheo))
    fichier.close()

    #ecriture du fichier pour la courbe Utau trio
    nomFicCourbeUtautrio = '../CourbeUtautrio_'+scheme+'.dat'
    fichier = open(nomFicCourbeUtautrio, 'a')
    fichier.write('%18.8f %18.8f\n' % (Re,Utrio))
    fichier.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    #recuperation des donnees
    scheme = args[0]
    U = args[1]
    mu,rho,h,L = properties()
    reh = Re(mu,rho,h,U)
    y = yPlus(scheme,U)
    Uth = Utheo(reh,mu,rho,h)
    Utr = Utrio(scheme,U)
    Uerror = Uerr(Uth,Utr)
    Pth = Ptheo(Uth,L,h)
    Ptr = Ptrio(scheme,U,L)
    Perror = Perr(Pth,Ptr)

    #ecriture des fichiers gnuplot et tableaux
    ecritureFichier(scheme,U,reh,y,Uth,Utr,Uerror,Pth,Ptr,Perror)
