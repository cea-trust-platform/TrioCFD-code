#Script de recuperation des donnees fraction massique et densite du melange

import optparse
import math

def analytic():
    nomFic = 'analyticResults.dat'

    fic = open(nomFic,'r')

    fichier = fic.readlines()

    ligne = fichier [0]
    tLigne = ligne.split()
    AnaMF1=float(tLigne[0])

    fic.close()
    return AnaMF1


def results():

    nomFicVDFMF1 = 'NoFlow-lami-espece_VDF_FRACTION_MASSIQUE0.son'
    FicVDFMF1 = open(nomFicVDFMF1,'r')
    fichier = FicVDFMF1.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VDFMF1=float(tLigne[1])
    FicVDFMF1.close()

    nomFicVEFMF1 = 'NoFlow-lami-espece_VEF_FRACTION_MASSIQUE0.son'
    FicVEFMF1 = open(nomFicVEFMF1,'r')
    fichier = FicVEFMF1.readlines()
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]
    ligne = fichier [-1]
    tLigne = ligne.split()
    VEFMF1=float(tLigne[1])
    FicVEFMF1.close()

    return VDFMF1,VEFMF1


def ecritureFichier(VDFMF1,VEFMF1,AnaMF1,ErrVDF,ErrVEF):
    #ecriture du fichier pour gnuplot
    nomFic = 'ligneTableau-specie1.dat'

    fichier = open(nomFic, 'w')
    fichier.write(' %8.4f %8.4f %8.4f %8.4f %8.4f\n' % (VDFMF1,VEFMF1,AnaMF1,ErrVDF,ErrVEF))
    fichier.close()


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    AnaMF1 = analytic()
    VDFMF1,VEFMF1 = results()
    ErrVDF=(abs(VDFMF1-AnaMF1)/AnaMF1)*100
    ErrVEF=(abs(VEFMF1-AnaMF1)/AnaMF1)*100

    #ecriture du fichier gnuplot
    ecritureFichier(VDFMF1,VEFMF1,AnaMF1,ErrVDF,ErrVEF)
