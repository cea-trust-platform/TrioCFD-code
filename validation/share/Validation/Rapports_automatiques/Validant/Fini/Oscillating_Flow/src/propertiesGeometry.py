#Script de recuperation des proprietes physiques et de la geometrie

import os, sys, math

def readPropertiesData(nomFic):
    #initialisation
    properties = {}
    properties['mu'] = -1
    properties['lambda'] = -1
    properties['cp'] = -1
    properties['rho'] = -1
    properties['beta_th'] = -1
    properties['alpha'] = -1
    properties['g'] = -1
    # ouverture des fichiers
    fic = open(nomFic,'r')

    #chaines = ["mu champ_uniforme",
    #       "lambda champ_uniforme",
    #       "beta_th champ_uniforme",
    #       "cp champ_uniforme",
    #       "rho champ_uniforme",
    #       "beta_th champ_uniforme",
    #       "vitesse champ_fonc_xyz",
    #       "paroi paroi_flux_impose champ_front_uniforme",
    #       "sonde_tpn ",
    #       "lire gravite"] # Texte a rechercher

    for ligne in fic:
        ligne = ligne.strip().lower()
        tLigne = ligne.split()
        if ligne.find('champ_uniforme')>-1:
            if ligne.startswith('mu'):
                properties['mu'] = float(tLigne[-1])
            elif ligne.startswith('lambda'):
                properties['lambda'] = float(tLigne[-1])
            elif ligne.startswith('beta_th'):
                properties['beta_th'] = float(tLigne[-1])
            elif ligne.startswith('cp'):
                properties['cp'] = float(tLigne[-1])
            elif ligne.startswith('rho'):
                properties['rho'] = float(tLigne[-1])
            #else:
            #       print 'Ligne avec champ_uniforme non reconnue : %s' % (ligne)
        elif ligne.find('gravite')>-1:
            if ligne.startswith('lire') or ligne.startswith('read'):
                dim = int(tLigne[2])
                g = 0.
                for i in range(dim):
                    gi = float(tLigne[i+3])
                    g += gi*gi
                properties['g'] = math.sqrt(g)
            #else:
            #       print 'Ligne avec gravite non reconnue : %s' % (ligne)


        #for chaine in chaines:
        #       if chaine.lower() in ligne:
        #               tLigne = ligne.split()
        #               elif chaine=="vitesse champ_fonc_xyz":
        #                       formule = tLigne[-1]
        #                       tFormule = formule.split('*')
        #                       U=float(tFormule[0])/2.000
        #               elif chaine=="paroi   paroi_flux_impose Champ_Front_Uniforme":
        #                       Qw=float(tLigne[-1])
        #               elif chaine=="sonde_tpn ":
        #                       Z=float(tLigne[-1])
        #                       D=2.*float(tLigne[9])
        #               else:
        #                       print 'chaine (%s) non reconnue' % (ligne)
        #       else:
        #               print 'chaine (%s) non reconnue' % (ligne)
    fic.close()
    properties['alpha'] = properties['lambda'] / (properties['rho']*properties['cp'])
    return properties


def ecritureFichier(properties):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f\n' % ( properties['mu'], properties['lambda'], properties['cp'], properties['rho'], properties['beta_th'], properties['alpha'], properties['g']))
    fichier.close()

def getPropertiesFromdat():
    properties = {}
    nomFichier = 'propertiesGeometry.dat'
    if os.path.isfile(nomFichier):
        #recupere les donnees du fichier
        f = open(nomFichier, 'r')
        lignes = f.readlines()
        f.close()
    else:
        print 'Erreur getPropertiesFromdat : fichier %s non trouve !' % (nomFichier)
        sys.exit()
    ligne = (lignes[0]).strip()
    tabLigne = ligne.split()
    ind = 0
    try:
        properties['mu'] = float(tabLigne[ind])
        ind += 1
        properties['lambda'] = float(tabLigne[ind])
        ind += 1
        properties['cp'] = float(tabLigne[ind])
        ind += 1
        properties['rho'] = float(tabLigne[ind])
        ind += 1
        properties['beta_th'] = float(tabLigne[ind])
        ind += 1
        properties['alpha'] = float(tabLigne[ind])
        ind += 1
        properties['g'] = float(tabLigne[ind])
    except IndexError:
        print 'Erreur getPropertiesFromdat : lecture element %d pour 0-%d elements...' % (ind, len(tabLigne)-1)
        sys.exit()
    except ValueError:
        print 'Erreur getPropertiesFromdat : lecture element %d n\'est pas un float (%s)...' % (ind, tabLigne[ind])
        sys.exit()
    return properties


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
        properties = readPropertiesData(nomFic)

        #ecriture du fichier
        ecritureFichier(properties)
    else:
        print 'Erreur propertiesGeometry : pas de fichier data trouve !'
