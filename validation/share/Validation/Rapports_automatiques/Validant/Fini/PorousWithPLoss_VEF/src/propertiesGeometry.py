#Script de recuperation des proprietes physiques et de la geometrie

import os, sys, math

def readPropertiesData(nomFic):
    #initialisation
    properties = {}
    properties['mu'] = -1
    properties['Cp'] = -1
    properties['rho'] = -1
    properties['beta'] = -1
    properties['lambda'] = -1
    properties['diam_hydr'] = -1
    properties['lambda_p'] = -1
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
        if ligne.startswith('mu'):
            properties['mu'] = float(tLigne[-1])
        elif ligne.startswith('beta_th'):
            properties['beta'] = float(tLigne[-1])
        elif ligne.startswith('lambda'):
            properties['lambda'] = float(tLigne[-1])
        elif ligne.startswith('cp'):
            properties['Cp'] = float(tLigne[-1])
        elif ligne.startswith('rho'):
            properties['rho'] = float(tLigne[-1])
        elif ligne.startswith('sources {'):
            properties['diam_hydr'] = float(tLigne[7])
            properties['lambda_p'] = float(tLigne[9])

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
    return properties


def ecritureFichier(properties):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f \n' % ( properties['mu'], properties['Cp'], properties['rho'], properties['beta'], properties['lambda'], properties['diam_hydr'], properties['lambda_p'] ))
    fichier.close()

if __name__ == '__main__':

    #recuperation du fichier data
    import glob
    listFics = glob.glob('*.data')
    if len(listFics)>0:
        nomFic = listFics[0]
        properties = readPropertiesData(nomFic)

        #ecriture du fichier
        ecritureFichier(properties)
    else:
        print 'Erreur propertiesGeometry : pas de fichier data trouve !'
