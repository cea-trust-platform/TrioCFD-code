#Script de construction des tableaux Euler_Explicit

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

def Re(mu,h,U,rho):
    reh = float(rho)*float(U)*2*float(h)/float(mu)
    return reh

def Utheo(reh,U):
    utau = float(U)*math.sqrt(12/reh)
    return utau

def Utrio():
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
    utau=float(tLigne[-1])
    ficUstar.close()

    return utau

def Uerr(Utheo,Utrio):
    Uerror = float(abs(Utheo-Utrio)/Utheo)*100
    return Uerror

def Ptheo(Utheo,L,h):
    Pth = float(Utheo*Utheo*L/(h/2))
    return Pth

def Ptrio(L):
    # ouverture des fichiers
    nomFic = os.popen('ls -rt *Pressure_Gradient_pb1_*').readlines()[-1]
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

def ecritureFichier(discr,scheme,U,Re,Utheo,Utrio,UError,Ptheo,Ptrio,PError):
    #ecriture du fichier pour la ligne du tableau
    nomFicLigne = 'ligneTableau.dat'
    fichier = open(nomFicLigne, 'w')
    fichier.write('%s %18.2f %e %e %18.2f %e %e %18.2f\n' % (U,Re,Utheo,Utrio,UError,Ptheo,Ptrio,PError))
    fichier.close()

    #ecriture du fichier pour la courbe dPtheo
    nomFicCourbedPtheo = 'CourbedPtheo.dat'
    fichier = open(nomFicCourbedPtheo, 'w')
    fichier.write('%18.8f %18.8f\n' % (Re,Ptheo))
    fichier.close()

    #ecriture du fichier pour la courbe dPtrio
    nomFicCourbedPtrio = 'CourbedPtrio.dat'
    fichier = open(nomFicCourbedPtrio, 'w')
    fichier.write('%18.8f %18.8f\n' % (Re,Ptrio))
    fichier.close()

    #ecriture du fichier pour la courbe Utau theo
    nomFicCourbeUtautheo = 'CourbeUtautheo.dat'
    fichier = open(nomFicCourbeUtautheo, 'w')
    fichier.write('%18.8f %18.8f\n' % (Re,Utheo))
    fichier.close()

    #ecriture du fichier pour la courbe Utau trio
    nomFicCourbeUtautrio = 'CourbeUtautrio.dat'
    fichier = open(nomFicCourbeUtautrio, 'w')
    fichier.write('%18.8f %18.8f\n' % (Re,Utrio))
    fichier.close()

def copieFichierMoyennesSpatiales(nomFic):
    #ecriture du fichier
    ficRead = open(nomFic,'r')
    ficWrite = open('Velocity_profile.dat', 'w')

    fin=False
    i=0
    while not fin:
        ligne = ficRead.readline()
        if not ligne:
            fin=True
        else:
            if ligne[0]!="#" and ligne[0]!=" " and ligne[0]!="":
                ficWrite.write('%s' % (ligne))

    ficWrite.close()
    ficRead.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    #recuperation des donnees
    discr = args[0]
    scheme = args[1]
    U = args[2]
    print(discr+' '+scheme+' '+U)
    mu,rho,h,L = properties()
    reh = Re(mu,h,U,rho)
    Uth = Utheo(reh,U)
    Utr = Utrio()
    Uerror = Uerr(Uth,Utr)
    Pth = Ptheo(Uth,L,h)
    Ptr = Ptrio(L)
    Perror = Perr(Pth,Ptr)

    #ecriture des fichiers gnuplot et tableaux
    ecritureFichier(discr,scheme,U,reh,Uth,Utr,Uerror,Pth,Ptr,Perror)

    #derniere ligne du ls
    ficLS = os.popen('ls -rt Moyennes_spatiales_vitesse_rho_mu_*')
    lignes = ficLS.readlines()
    derLigne = lignes[-1]
    #suppression du \n en fin de nom
    nomFic = derLigne[:len(derLigne)-1]

    #copie du fichier MoyennesSpatiales* dans un fichier generique VelocityProfile
    copieFichierMoyennesSpatiales(nomFic)
