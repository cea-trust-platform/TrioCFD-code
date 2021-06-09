#Script de recuperation des donnees fraction massique et densite du melange

import optparse
import math

def analytic():
    nomFic = 'analyticResults.dat'

    fic = open(nomFic,'r')

    fichier = fic.readlines()

    ligne = fichier [0]
    tLigne = ligne.split()
    AnaMF2=float(tLigne[1])

    fic.close()
    return AnaMF2


def results():

    nomFicVDFMF2 = 'NoFlow-lami-espece_VDF_FRACTION_MASSIQUE1.son'
    FicVDFMF2 = open(nomFicVDFMF2,'r')
    fichier = FicVDFMF2.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VDFMF2=float(tLigne[1])
    FicVDFMF2.close()

    nomFicVEFMF2 = 'NoFlow-lami-espece_VEF_FRACTION_MASSIQUE1.son'
    FicVEFMF2 = open(nomFicVEFMF2,'r')
    fichier = FicVEFMF2.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VEFMF2=float(tLigne[1])
    FicVEFMF2.close()

    return VDFMF2,VEFMF2


def ecritureFichier(VDFMF2,VEFMF2,AnaMF2,ErrVDF,ErrVEF):
    #ecriture du fichier pour gnuplot
    nomFic = 'ligneTableau-specie2.dat'

    fichier = open(nomFic, 'w')
    fichier.write(' %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (VDFMF2,VEFMF2,AnaMF2,ErrVDF,ErrVEF))
    fichier.close()


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    AnaMF2 = analytic()
    VDFMF2,VEFMF2 = results()
    ErrVDF=(abs(VDFMF2-AnaMF2)/AnaMF2)*100
    ErrVEF=(abs(VEFMF2-AnaMF2)/AnaMF2)*100

    #ecriture du fichier gnuplot
    ecritureFichier(VDFMF2,VEFMF2,AnaMF2,ErrVDF,ErrVEF)
