#Script de recuperation des donnees fraction massique et densite du melange

import optparse
import math

def analytic():
    nomFic = 'analyticResults.dat'

    fic = open(nomFic,'r')

    fichier = fic.readlines()

    ligne = fichier [0]
    tLigne = ligne.split()
    AnaDensity=float(tLigne[2])

    fic.close()
    return AnaDensity


def results():

    nomFicVDFDensity1 = 'NoFlow-lami-espece_VDF_MASSE_VOLUMIQUE1.son'
    FicVDFDensity1 = open(nomFicVDFDensity1,'r')
    fichier = FicVDFDensity1.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VDFDensity1=float(tLigne[1])
    FicVDFDensity1.close()

    nomFicVDFDensity2 = 'NoFlow-lami-espece_VDF_MASSE_VOLUMIQUE2.son'
    FicVDFDensity2 = open(nomFicVDFDensity2,'r')
    fichier = FicVDFDensity2.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VDFDensity2=float(tLigne[1])
    FicVDFDensity2.close()

    nomFicVEFDensity1 = 'NoFlow-lami-espece_VEF_MASSE_VOLUMIQUE1.son'
    FicVEFDensity1 = open(nomFicVEFDensity1,'r')
    fichier = FicVEFDensity1.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VEFDensity1=float(tLigne[1])
    FicVEFDensity1.close()

    nomFicVEFDensity2 = 'NoFlow-lami-espece_VEF_MASSE_VOLUMIQUE2.son'
    FicVEFDensity2 = open(nomFicVEFDensity2,'r')
    fichier = FicVEFDensity2.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VEFDensity2=float(tLigne[1])
    FicVEFDensity2.close()

    return VDFDensity1,VDFDensity2,VEFDensity1,VEFDensity2


def ecritureFichier(VDFDensity,VEFDensity,AnaDensity,ErrVDF,ErrVEF):
    #ecriture du fichier pour gnuplot
    nomFic = 'ligneTableau-density.dat'

    fichier = open(nomFic, 'w')
    fichier.write(' %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (VDFDensity,VEFDensity,AnaDensity,ErrVDF,ErrVEF))
    fichier.close()


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    AnaDensity = analytic()
    VDFDensity1,VDFDensity2,VEFDensity1,VEFDensity2 = results()
    VDFDensity=(VDFDensity1+VDFDensity2)/2
    VEFDensity=(VEFDensity1+VEFDensity2)/2
    ErrVDF=(abs(VDFDensity-AnaDensity)/AnaDensity)*100
    ErrVEF=(abs(VEFDensity-AnaDensity)/AnaDensity)*100

    #ecriture du fichier gnuplot
    ecritureFichier(VDFDensity,VEFDensity,AnaDensity,ErrVDF,ErrVEF)
