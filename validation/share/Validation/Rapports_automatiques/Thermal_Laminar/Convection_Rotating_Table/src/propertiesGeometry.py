#Script de recuperation des propretes physiques et de la geometrie

import os, sys

def readPropertiesData(nomFic):
    #initialisation
    D = -1
    mu = -1
    lamb = -1
    Cp = -1
    rho = -1
    beta_th = -1
    alpha = -1
    U = -1
    Qw = -1
    Z = -1
    beta_th = -1
    # ouverture des fichiers
    fic = open(nomFic,'r')

    chaines = ["mu Champ_Uniforme",
               "lambda Champ_Uniforme",
               "beta_th Champ_Uniforme",
               "Cp Champ_Uniforme",
               "rho Champ_Uniforme",
               "beta_th Champ_Uniforme",
               "vitesse champ_fonc_xyz",
               "paroi   paroi_flux_impose Champ_Front_Uniforme",
               "sonde_tpn "] # Texte a rechercher

    for ligne in fic:
        for chaine in chaines:
            if chaine in ligne:
                tLigne = ligne.split()
                if chaine=="mu Champ_Uniforme":
                    mu=float(tLigne[-1])
                elif chaine=="lambda Champ_Uniforme":
                    lamb=float(tLigne[-1])
                elif chaine=="Cp Champ_Uniforme":
                    Cp=float(tLigne[-1])
                elif chaine=="rho Champ_Uniforme":
                    rho=float(tLigne[-1])
                elif chaine=="vitesse champ_fonc_xyz":
                    formule = tLigne[-1]
                    tFormule = formule.split('*')
                    U=float(tFormule[0])/2.000
                elif chaine=="paroi   paroi_flux_impose Champ_Front_Uniforme":
                    Qw=float(tLigne[-1])
                elif chaine=="sonde_tpn ":
                    Z=float(tLigne[-1])
                    D=2.*float(tLigne[9])
                elif chaine=="beta_th Champ_Uniforme":
                    beta_th =float(tLigne[-1])
    fic.close()
    alpha = lamb / (rho*Cp)
    return D,mu,lamb,Cp,rho,beta_th,alpha,U,Qw,Z


def ecritureFichier(D,mu,lamb,Cp,rho,beta_th,alpha,U,Qw,Z):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.19f %18.19f %18.19f %18.19f %18.19f %18.19f %18.19f %18.19f %18.19f %18.19f\n' % (D,mu,lamb,Cp,rho,beta_th,alpha,U,Qw,Z))
    fichier.close()

def getPropertiesFromdat():
    nomFichier = 'propertiesGeometry.dat'
    if os.path.isfile(nomFichier):
            #recupere les donnees du fichier
        f = open(nomFichier, 'r')
        lignes = f.readlines()
        f.close()
    else:
        print('Erreur getPropertiesFromdat : fichier %s non trouve !' % (nomFichier))
        sys.exit()
    ligne = (lignes[0]).strip()
    tabLigne = ligne.split()
    ind = 0
    try:
        D       = float(tabLigne[ind])
        ind += 1
        mu      = float(tabLigne[ind])
        ind += 1
        lamb    = float(tabLigne[ind])
        ind += 1
        Cp      = float(tabLigne[ind])
        ind += 1
        rho     = float(tabLigne[ind])
        ind += 1
        beta_th = float(tabLigne[ind])
        ind += 1
        alpha   = float(tabLigne[ind])
        ind += 1
        U       = float(tabLigne[ind])
        ind += 1
        Qw      = float(tabLigne[ind])
        ind += 1
        Z       = float(tabLigne[ind])
    except IndexError:
        print('Erreur getPropertiesFromdat : lecture element %d pour 0-%d elements...' % (ind, len(tabLigne)-1))
        sys.exit()
    except ValueError:
        print('Erreur getPropertiesFromdat : lecture element %d n\'est pas un float (%s)...' % (ind, tabLigne[ind]))
        sys.exit()
    return D,mu,lamb,Cp,rho,beta_th,alpha,U,Qw,Z


def nombreDePrandtl(mu, lamb, Cp):
    Pr = mu / (lamb / Cp)
    return Pr

def nombreDeRayleigh(g, beta_th, rho, mu, alpha, Tparoi, Tfluide, Lcarac):
    Ra = g*beta_th * (Tparoi - Tfluide) * Lcarac * Lcarac * Lcarac /(mu/rho * alpha)
    return Ra

if __name__ == '__main__':

    #recuperation du fichier data
    import glob
    #derniere ligne du ls
    #ficLS = os.popen('ls *.data')
    #lignes = ficLS.readlines()
    #Ligne = lignes[0]
    #suppression du \n en fin de nom
    #nomFic = Ligne[:len(Ligne)-1]
    listFics = glob.glob('*.data')
    if len(listFics)>0:
        nomFic = listFics[0]
        D,mu,lamb,Cp,rho,beta_th,alpha,U,Qw,Z = readPropertiesData(nomFic)

        #ecriture du fichier
        ecritureFichier(D,mu,lamb,Cp,rho,beta_th,alpha,U,Qw,Z)
    else:
        print('Erreur propertiesGeometry : pas de fichier data trouve !')
