#Script de recuperation des proprietes physiques et de la geometrie

import os, sys, math

def readPropertiesData(nomFic):
    #initialisation
    properties = {}
    properties['mu'] = -1
    properties['rho_f'] = -1
    properties['lambda_f'] = -1
    properties['Cp'] = -1
    properties['beta_th'] = -1
    properties['rho_a'] = -1
    properties['lambda_a'] = -1
    # ouverture des fichiers
    fic = open(nomFic,'r')

#

    for ligne in fic:
        ligne = ligne.strip().lower()
        tLigne = ligne.split()
        if ligne.find('champ_uniforme')>-1:
            if ligne.startswith('mu'):
                properties['mu'] = float(tLigne[-1])
            elif ligne.startswith('rho champ_uniforme 1 27'):
                properties['rho_f'] = float(tLigne[-1])
            elif ligne.startswith('lambda champ_uniforme 1 0'):
                properties['lambda_f'] = float(tLigne[-1])
            elif ligne.startswith('cp'):
                properties['Cp'] = float(tLigne[-1])
            elif ligne.startswith('beta_th'):
                properties['beta_th'] = float(tLigne[-1])
            elif ligne.startswith('rho champ_uniforme 1 1'):
                properties['rho_a'] = float(tLigne[-1])
            elif ligne.startswith('lambda champ_uniforme 1 15'):
                properties['lambda_a'] = float(tLigne[-1])

    fic.close()
    return properties


def ecritureFichier(properties):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')
    fichier.write('%18.7f %18.3f %18.3f %18.3f %18.7f %18.3f %18.3f\n' % ( properties['mu'], properties['rho_f'], properties['lambda_f'], properties['Cp'], properties['beta_th'],  properties['rho_a'], properties['lambda_a'] ))
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
        print('Erreur getPropertiesFromdat : fichier %s non trouve !' % (nomFichier))
        sys.exit()
    ligne = (lignes[0]).strip()
    tabLigne = ligne.split()
    ind = 0
    try:
        properties['mu'] = float(tabLigne[ind])
        ind += 1
        properties['rho_f'] = float(tabLigne[ind])
        ind += 1
        properties['lambda_f'] = float(tabLigne[ind])
        ind += 1
        properties['Cp'] = float(tabLigne[ind])
        ind += 1
        properties['beta_th'] = float(tabLigne[ind])
        ind += 1
        properties['rho_a'] = float(tabLigne[ind])
        ind += 1
        properties['lambda_a'] = float(tabLigne[ind])
    except IndexError:
        print('Erreur getPropertiesFromdat : lecture element %d pour 0-%d elements...' % (ind, len(tabLigne)-1))
        sys.exit()
    except ValueError:
        print('Erreur getPropertiesFromdat : lecture element %d n\'est pas un float (%s)...' % (ind, tabLigne[ind]))
        sys.exit()
    return properties



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
        print('Erreur propertiesGeometry : pas de fichier data trouve !')
