import os,os.path, sys, math
############
############    Calcul de la des valeurs moyenne sur un segment (selon x)
############
############
## -- somme des valeurs sur toutes les faces du segment (selon x)
## -- sonde_A.son ==> Som(rho * T * U)
## -- sonde_B.son ==> Som(rho * U)
## -- sonde_C.son ==> Som(rho * T)
##
## - NN ==> nombre de faces du maillage selon x
##
##      Pmoyen = (1/NN) * Som(P)
##
##      Tmoyen = sonde_A / sonde_B
##
##      Rhomoyen = (1/NN) * sonde_B * sonde_C / sonde_A
##
##      Vmoyen = sonde_A / sonde_C
##
###

yy = [0., 0.2, 0.4, 0.6, 0.87406]
listeCas = ['0', '02', '04', '06', '087']

def calculerSomme(nomFichier, offset, isVecteur):
    f = open(nomFichier, 'r')
    # on saute les 5 premieres lignes
    for i in range(6):
        ligne = f.readline()
    ligneTab = ligne.strip().split()
    #calcule de la somme des valeurs de la ligne (sauf 1e valeur = temps)
    if isVecteur:
        #si c'est un vecteur : on ne recupere qu'un valeur sur 3
        ligneTab = ligneTab[offset::3]
    else:
        ligneTab = ligneTab[offset:]
    somme = 0.
    for val in ligneTab:
        somme += float(val)
    nbElem = len(ligneTab)
    return somme, nbElem

def calculerPmoy():
    fic = open('pmoyen.dat','w')
    for i, cas in enumerate(listeCas):
        nomFic = 'test_PMOY%s.son' % (cas)
        somme,nbElem = calculerSomme(nomFic, 1, False)
        moy = somme/nbElem
        fic.write('%18.2f %18.2f\n' % (yy[i], abs(moy)))
    fic.close()

def calculerSonde(type):
    vals = []
    for cas in ['0', '02', '04', '06', '087']:
        nomFic = 'test_SONDE_%s%s.son' % (type,cas)
        if type=='C':
            somme,nbElem = calculerSomme(nomFic, 1, False)
        else:
            somme,nbElem = calculerSomme(nomFic, 2, True)
        vals.append(somme)
    return vals,nbElem


def calculerMoySondes(A, B, C, nbElem,TBU):
    Tmoy = 0.
    Vmoy = 0.
    Rmoy = 0.
    fa = open('averages.dat','w')
    fb = open('err.dat','w')

#
    fichier = TBU.readlines()
    l = fichier [0]
    tLi = l.split()
    ttb = float(tLi[0])
#
    for p in range(len(listeCas)):
        Tmoy = A[p] / B[p]
        Vmoy = A[p] / C[p]
        Rmoy = B[p] * C[p] / ( nbElem * A[p] )
        y = yy[p]
        if p!=5:
            tt = Tmoy
            err = abs(100*(ttb-tt))/max(ttb,tt)
            fb.write('%18.2f\n' % (err))
        fa.write('%18.2f %18.2f %18.2f %18.2f\n' % (y, Tmoy, Vmoy, Rmoy))
    fa.close()

def calculerSondes(TBU):
    A,nbElem = calculerSonde('A')
    B,nbElem = calculerSonde('B')
    C,nbElem = calculerSonde('C')
    calculerMoySondes(A, B, C, nbElem,TBU)

if __name__=='__main__':
    args = sys.argv
    TB = args[21]
    TBU = open(TB,'r')
    #calcul pour PMoy
    calculerPmoy()
    #calcul pour les sondes 1, B, C
    calculerSondes(TBU)
