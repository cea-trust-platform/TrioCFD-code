import os,os.path, sys, math

def getNormePression(nomFichier, numPoint, fromTime):
    temps = []
    pression = []
    if os.path.isfile(nomFichier):
        f = open(nomFichier, 'r')
        lignes = f.readlines()
        f.close()
        for ligne in lignes:
            ligne = ligne.strip()
            if ligne[0]!='#':
                tLigne = ligne.split()
                valtemps = float(tLigne[0])
                valpression = float(tLigne[numPoint])
                if valtemps>=fromTime:
                    temps.append(valtemps)
                    pression.append(valpression)
    else:
        print 'Erreur : fichier non trouve %s' % (nomFichier)
    return temps, pression

def calculerMediane(tabOrig):
    #duplication du tableau d'origine
    tab = []
    for el in tabOrig:
        tab.append(el)
    mediane = 0
    n = len(tab)
    i = 0
    while i<n:
        j = 0
        while j<n-1-i:
            if tab[j]>tab[j+1]:
                tmp = tab[j+1]
                tab[j+1] = tab[j]
                tab[j] = tmp
            j += 1
        i += 1
    mediane = tab[(n+1)/2]
    return mediane

def getPeriodes(temps, tab, mediane):
    periodes = []
    nbCoupes = -1
    i = 0
    n = len(tab)
    periodeMoy = 0.
    while i<n-1:
        etati = (tab[i]<mediane)
        etatip1 = (tab[i+1]<mediane)
        if etati!=etatip1:
            #on vient de traverser la mediane
            nbCoupes += 1
            #print 'coupe mediane a temps= %f (%d)' % (temps[i], nbCoupes)
            if nbCoupes==0:
                temps0 = temps[i]
            else:
                if (nbCoupes % 2)==0:
                    periode = temps[i] - temps0
                    periodes.append(periode)
                    temps0 = temps[i]
                    periodeMoy += periode
        i += 1
    periodeMoy = periodeMoy / (len(periodes))
    return periodes, periodeMoy

if __name__=='__main__':
    #recuperation des parametres passes en ligne de commande
    args = sys.argv
    if len(args)!=4:
        print 'Erreur sur le nb d\'arguments fournis : Usage\npython calculerFrequence nomFichier.sonde numPoint fromTime'
        sys.exit()

    nomFic = args[1]
    numPoint = int(args[2])
    fromTime = float(args[3])
    temps, pression = getNormePression(nomFic,numPoint, fromTime)
    mediane = calculerMediane(pression)
    #print 'mediane = %f' % (mediane)

    periodes, periodeMoy = getPeriodes(temps, pression, mediane)
    #print 'Periodes=%s, moyenne=%f' % (periodes, periodeMoy)

    frequence = 1./periodeMoy
    print '%f' % (frequence)
