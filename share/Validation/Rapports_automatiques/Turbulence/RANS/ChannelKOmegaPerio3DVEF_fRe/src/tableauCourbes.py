# Script de construction des tableaux Euler_Explicit et Crank_Nicholson et leurs courbes

import optparse
import math
import os
from glob import glob


# Fonction de recuperation de la colonne y+
def yPlus(scheme, U):

    # ouverture des fichiers
    nomFicUstar = 'canal_turbu_'+scheme+'_pb1_Ustar.face'

    with open(nomFicUstar, 'r') as ficUstar:

        # lecture de ligne -> entetes
        fichier = ficUstar.readlines()

        # tant que la derniere ligne est vide
        while fichier[-1] == "" or fichier[-1] == "\n":
            del fichier[-1]

        ligne = fichier[-1]
        tLigne = ligne.split()
        y = float(tLigne[3])

    return y


def properties():
    # ouverture des fichiers
    nomFic = 'propertiesGeometry.dat'

    with open(nomFic, 'r') as fic:

        # lecture de ligne -> entetes
        fichier = fic.readlines()

        ligne = fichier[0]
        tLigne = ligne.split()
        mu = float(tLigne[0])
        rho = float(tLigne[1])
        h = float(tLigne[2])
        L = float(tLigne[3])

    return mu, rho, h, L


def Re(mu, rho, h, U):
    reh = float(rho)*float(U)*2*float(h)/float(mu)
    return reh


def Utheo(reh, mu, rho, h):
    rebulk = float(reh/4)
    retau = 0.1754*math.pow(float(rebulk), (7.0/8.0))
    utau = mu*retau/rho/(h/2)
    return utau


def Utrio(scheme, U):
    # ouverture des fichiers
    nomFicUstar = 'canal_turbu_'+scheme+'_pb1_Ustar.face'

    with open(nomFicUstar, 'r') as ficUstar:

        # lecture de ligne -> entetes
        fichier = ficUstar.readlines()

        # tant que la derniere ligne est vide
        while fichier[-1] == "" or fichier[-1] == "\n":
            del fichier[-1]

        ligne = fichier[-1]
        tLigne = ligne.split()
        utau = float(tLigne[5])

    return utau


def Uerr(Utheo, Utrio):
    Uerror = float(abs(Utheo-Utrio)/Utheo)*100
    return Uerror


def Ptheo(Utheo, L, h):
    Pth = float(Utheo*Utheo*L/(h/2))
    return Pth


def Ptrio(scheme, U, L):
    # ouverture des fichiers
    nomFic = os.popen('ls -rt *Pressure_Gradient_pb1_Entree').readlines()[-1]
    nomFic = nomFic[:len(nomFic)-1]  # Suppress /n

    with open(nomFic, 'r') as fic:

        # lecture de ligne -> entetes
        fichier = fic.readlines()

        # tant que la derniere ligne est vide
        while fichier[-1] == "" or fichier[-1] == "\n":
            del fichier[-1]

        ligne = fichier[-1]
        tLigne = ligne.split()
        ptrio = float(tLigne[1])*float(L)

    return ptrio


def Perr(Ptheo, Ptrio):
    Uerror = float(abs(Ptheo-Ptrio)/Ptheo)*100
    return Uerror


def ecritureFichier(scheme, U, Re, y, Utheo, Utrio, UError, Ptheo, Ptrio, PError):
    # ecriture du fichier pour la ligne du tableau
    nomFicLigne = 'ligneTableau.dat'
    with open(nomFicLigne, 'w') as fichier:
        fichier.write('%s %18.2f %18.2f %e %e %18.2f %e %e %18.2f\n' % (U, Re, y, Utheo, Utrio, UError, Ptheo, Ptrio, PError))

    # ecriture du fichier pour la courbe dPtheo
    nomFicCourbedPtheo = '../CourbedPtheo_'+scheme+'.dat'
    with open(nomFicCourbedPtheo, 'a') as fichier:
        fichier.write('%18.8f %18.8f\n' % (Re, Ptheo))

    # ecriture du fichier pour la courbe dPtrio
    nomFicCourbedPtrio = '../CourbedPtrio_'+scheme+'.dat'
    with open(nomFicCourbedPtrio, 'a') as fichier:
        fichier.write('%18.8f %18.8f\n' % (Re, Ptrio))

    # ecriture du fichier pour la courbe Utau theo
    nomFicCourbeUtautheo = '../CourbeUtautheo_'+scheme+'.dat'
    with open(nomFicCourbeUtautheo, 'a') as fichier:
        fichier.write('%18.8f %18.8f\n' % (Re, Utheo))

    # ecriture du fichier pour la courbe Utau trio
    nomFicCourbeUtautrio = '../CourbeUtautrio_'+scheme+'.dat'
    with open(nomFicCourbeUtautrio, 'a') as fichier:
        fichier.write('%18.8f %18.8f\n' % (Re, Utrio))


def copieFichierMoyennesSpatiales(nomFic):
    # ecriture du fichier de profil de vitesse
    with open(nomFic, 'r') as ficRead:
        with open('Velocity_profile.dat', 'w') as ficWrite:
            fin = False
            while not fin:
                ligne = ficRead.readline()
                if not ligne:
                    fin = True
                else:
                    if ligne[0] != "#" and ligne[0] != " " and ligne[0] != "":
                        ficWrite.write('%s' % (ligne))


def copieFichierMoyennesSpatiales_nut(nomFic):
    # ecriture du fichier de profil de viscosite turbulente
    with open(nomFic, 'r') as ficRead:
        with open('Turbulent_viscosity_profile.dat', 'w') as ficWrite:
            fin = False
            while not fin:
                ligne = ficRead.readline()
                if not ligne:
                    fin = True
                else:
                    if ligne[0] != "#" and ligne[0] != " " and ligne[0] != "":
                        ficWrite.write('%s' % (ligne))


def dimensionLess_velocity_profile(profil, Utr, nu):
    with open(profil, 'r') as ficRead:
        with open('Velocity_profile_dimensionLess.dat', 'w') as ficWrite:
            fin = False
            while not fin:
                ligne = ficRead.readline()
                if not ligne:
                    fin = True
                else:
                    tligne = ligne.split()
                    ficWrite.write('%18.8f %18.8f\n' % (float(tligne[0])*Utr/nu, float(tligne[1])/Utr))


def dimensionLess_turbulent_viscosity_profile(profil, Utr, nu):
    with open(profil, 'r') as ficRead:
        with open('Turbulent_viscosity_profile_dimensionLess.dat', 'w') as ficWrite:
            fin = False
            while not fin:
                ligne = ficRead.readline()
                if not ligne:
                    fin = True
                else:
                    tligne = ligne.split()
                    ficWrite.write('%18.8f %18.8f\n' % (float(tligne[0])*Utr/nu, float(tligne[1])/nu))


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    # recuperation des donnees
    scheme = args[0]
    U = args[1]
    mu, rho, h, L = properties()
    nu = float(mu/rho)
    reh = Re(mu, rho, h, U)
    y = yPlus(scheme, U)
    Uth = Utheo(reh, mu, rho, h)
    Utr = Utrio(scheme, U)
    Uerror = Uerr(Uth, Utr)
    Pth = Ptheo(Uth, L, h)
    Ptr = Ptrio(scheme, U, L)
    Perror = Perr(Pth, Ptr)

    # ecriture des fichiers gnuplot et tableaux
    ecritureFichier(scheme, U, reh, y, Uth, Utr, Uerror, Pth, Ptr, Perror)
    ###################################################################################

    # TRACE DES PROFILS DE VITESSE
    # derniere ligne du ls pour recuperer moyennes de vitesse.
    racine = "Moyennes_spatiales_vitesse_rho_mu_"
    maxfic = max([int(nn.replace(racine, "").split(".")[0]) for nn in glob(f"{racine}*")])
    nomFic = glob(f"{racine}{maxfic}*")[0]

    # copie du fichier MoyennesSpatiales* dans un fichier generique VelocityProfile
    copieFichierMoyennesSpatiales(nomFic)

    # Ecriture fichier profil vitesse adimensionne
    dimensionLess_velocity_profile('Velocity_profile.dat', Utr, nu)
    ###################################################################################

    # TRACE DES PROFILS DE VISCOSITE TURBULENTE
    # derniere ligne du ls pour recuperer moyennes de viscosite turbulente.
    racine = "Moyennes_spatiales_nut_"
    maxfic = max([int(nn.replace(racine, "").split(".")[0]) for nn in glob(f"{racine}*")])
    nomFic = glob(f"{racine}{maxfic}*")[0]

    # copie du fichier MoyennesSpatiales* dans un fichier generique VelocityProfile
    copieFichierMoyennesSpatiales_nut(nomFic)

    # Ecriture fichier profil viscosite turbulente adimensionnee
    dimensionLess_turbulent_viscosity_profile('Turbulent_viscosity_profile.dat', Utr, nu)
