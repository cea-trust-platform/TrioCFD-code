#Script de construction des fichiers de courbes

import optparse
import math
import os

def ecritureFichier(lignes):
    #courbe dPtheo
    i=0
    while i<len(lignes):
        ligne = lignes[i]
        #suppression du \n en fin de nom
        ligne = ligne[:len(ligne)-1]
        tLigne=ligne.split()
        nomDossier=tLigne[-1]
        elements=nomDossier.split('_')

        #ecriture du fichier pour la courbe dPtheo
        nomFicCourbedPtheo = 'CourbedPtheo_'+elements[0]+'_'+elements[1]+'.dat'
        fichier = open(nomFicCourbedPtheo, 'a')

        ficRead = open(nomDossier+'/CourbedPtheo.dat')
        lgn = ficRead.readline()

        fichier.write('%s' % (lgn))
        ficRead.close()
        fichier.close()
        i=i+1

    #courbe dPtrio
    i=0
    while i<len(lignes):
        ligne = lignes[i]
        #suppression du \n en fin de nom
        ligne = ligne[:len(ligne)-1]
        tLigne=ligne.split()
        nomDossier=tLigne[-1]
        elements=nomDossier.split('_')

        #ecriture du fichier pour la courbe dPtrio
        nomFicCourbedPtheo = 'CourbedPtrio_'+elements[0]+'_'+elements[1]+'.dat'
        fichier = open(nomFicCourbedPtheo, 'a')

        ficRead = open(nomDossier+'/CourbedPtrio.dat')
        lgn = ficRead.readline()

        fichier.write('%s' % (lgn))
        ficRead.close()
        fichier.close()
        i=i+1

    #courbe Utau theo
    i=0
    while i<len(lignes):
        ligne = lignes[i]
        #suppression du \n en fin de nom
        ligne = ligne[:len(ligne)-1]
        tLigne=ligne.split()
        nomDossier=tLigne[-1]
        elements=nomDossier.split('_')

        #ecriture du fichier pour la courbe Utau theo
        nomFicCourbedPtheo = 'CourbeUtautheo_'+elements[0]+'_'+elements[1]+'.dat'
        fichier = open(nomFicCourbedPtheo, 'a')

        ficRead = open(nomDossier+'/CourbeUtautheo.dat')
        lgn = ficRead.readline()

        fichier.write('%s' % (lgn))
        ficRead.close()
        fichier.close()
        i=i+1

    #courbe Utau trio
    i=0
    while i<len(lignes):
        ligne = lignes[i]
        #suppression du \n en fin de nom
        ligne = ligne[:len(ligne)-1]
        tLigne=ligne.split()
        nomDossier=tLigne[-1]
        elements=nomDossier.split('_')

        #ecriture du fichier pour la courbe Utau trio
        nomFicCourbedPtheo = 'CourbeUtautrio_'+elements[0]+'_'+elements[1]+'.dat'
        fichier = open(nomFicCourbedPtheo, 'a')

        ficRead = open(nomDossier+'/CourbeUtautrio.dat')
        lgn = ficRead.readline()

        fichier.write('%s' % (lgn))
        ficRead.close()
        fichier.close()
        i=i+1


if __name__ == '__main__':

    os.popen('rm -f Courbe*')
    ficLS = os.popen('ls -l | grep ^d')
    lignes = ficLS.readlines()

    #ecriture des fichiers gnuplot
    ecritureFichier(lignes)
