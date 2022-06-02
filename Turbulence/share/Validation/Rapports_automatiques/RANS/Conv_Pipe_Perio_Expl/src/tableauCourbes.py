#Script de construction des tableaux et des courbes
import os

import optparse
import math

def properties():
    # ouverture des fichiers
    nomFic = 'propertiesGeometry.dat'

    fic = open(nomFic,'r')

    # lecture de ligne -> entetes
    fichier = fic.readlines()

    ligne = fichier [0]
    tLigne = ligne.split()
    D=float(tLigne[0])
    mu=float(tLigne[1])
    lamb=float(tLigne[2])
    Cp=float(tLigne[3])
    rho=float(tLigne[4])
    U=float(tLigne[5])
    Qw=float(tLigne[6])
    Z=float(tLigne[7])

    fic.close()
    return D,mu,lamb,Cp,rho,U,Qw,Z

#fonction de recuperation de la colonne y+
def yPlus():
    # ouverture des fichiers
    nomFicUstar = 'Conduite_pb_Ustar.face'

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #en periodique on prend la valeur moyenne calculee en fin de fichier
    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    y=float(tLigne[3])
    ficUstar.close()

    return y

def uPlus():
    # ouverture des fichiers
    nomFicUstar = 'Conduite_pb_Ustar.face'

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    u=float(tLigne[1])
    ficUstar.close()

    return u

def Utrio():
    # ouverture des fichiers
    nomFicUstar = 'Conduite_pb_Ustar.face'

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    utau=float(tLigne[5])
    ficUstar.close()

    return utau

#fonction permettant de calculer la vitesse moyenne calculee par le code (differente de Ucible?)
def UbulkCalc(rho,D):
    #calcul de la la surface : attention dans RhoU_sortie surface double
    Surf=3.14159*float(D)*float(D)/4.
    #ouverture des fichiers
    nomFicUbulk = 'RhoU_sortie'
    ficUbulk = open(nomFicUbulk,'r')

    fichier = ficUbulk.readlines()
    #derniere ligne
    ligne = fichier[-1]
    tLigne = ligne.split()
    Ucalc=tLigne[1]
    Ucalc=float(Ucalc)/float(rho)/(2.*float(Surf))
    ficUbulk.close()
    return Ucalc

def Tparoi(Z):
    # recherche de Tparoi au meme endroit que la sonde pour tracer les profils
    # ouverture des fichiers
    nomFicTw = 'Conduite_pb_Nusselt.face'

    ficTw = open(nomFicTw,'r')

    tFic = ficTw.readlines()
    commentaire=""
    i=len(tFic)

    valz=Z
    deltaz=Z
    #determination de la valeur la plus proche de Z (par valeur inf)
    while commentaire!="----------------":
        ligne=tFic[i-1]
        if ligne != "\n" and ligne != "":
            tLigne = ligne.split("|")
            commentaire = tLigne[0]
            if commentaire!="----------------":
                #on eclate la ligne en un tableau de valeurs
                tLigne=ligne.split("|")
                # on lit la troisieme colonne (z)
                z=float(tLigne[2])
                if (Z-z) < abs(deltaz) and (Z-z)>0 :
                    valz = z
                    deltaz = Z-z
        i=i-1
    #on cherche la valeur x la plus proche de la paroi (X=-1) pour Z
    #attention on predefinit X=-1 (lire dans sonde, ce serait mieux)
    commentaire=""
    i=len(tFic)
    X=-1.
    valx=X
    deltaX=X
    while commentaire!="----------------":
        ligne=tFic[i-1]
        if ligne != "\n" and ligne != "":
            tLigne = ligne.split("|")
            commentaire = tLigne[0]
            if commentaire!="----------------":
                #on eclate la ligne en un tableau de valeurs
                tLigne=ligne.split("|")
                if float(tLigne[2])==valz :
                    x = float(tLigne[0])
                    if abs(X-x) <= abs(deltaX) :
                        valx = x
                        deltaX = X-x
        i=i-1

    #moyenne des valeurs Tparoi sur z = Z et x = X
    commentaire=""
    i=len(tFic)
    cumul=0
    nbTemps=0
    while commentaire!="----------------":
        ligne=tFic[i-1]
        if ligne != "\n" and ligne != "":
            tLigne = ligne.split("|")
            commentaire = tLigne[0]
            if commentaire!="----------------":
                #on eclate la ligne en un tableau de valeurs
                tLigne=ligne.split("|")
                if float(tLigne[2])==valz and float(tLigne[0])==valx:
                    # on lit la derniere valeur (Tp equiv)
                    T=float(tLigne[-1])
                    cumul = cumul + T
                    nbTemps = nbTemps + 1
        i=i-1

    Tw = cumul/float(nbTemps)

    ficTw.close()

    return Tw

def Ttotal():
    #en periodique on considere qu'on peut faire la moyenne sur toute la longueur
    # ouverture des fichiers
    nomFicTw = 'Conduite_pb_Nusselt.face'

    ficTw = open(nomFicTw,'r')

    tFic = ficTw.readlines()
    commentaire=""
    i=len(tFic)
    cumul=0
    nbTemps=0
    while commentaire!="----------------":
        ligne=tFic[i-1]
        if ligne != "\n" and ligne != "":
            tLigne = ligne.split("|")
            commentaire = tLigne[0]
            if commentaire!="----------------":
                #on eclate la ligne en un tableau de valeurs
                tLigne=ligne.split("|")
                # on lit la derniere valeur
                T=float(tLigne[-1])
                cumul = cumul + T
                nbTemps = nbTemps + 1
        i=i-1

    Tw = cumul/float(nbTemps)

    ficTw.close()

    return Tw

def Tbulk():
    #ouverture des fichiers
    nomFicTbulk = 'Tmoyen_sortie'

    ficTbulk = open(nomFicTbulk,'r')

    fichier = ficTbulk.readlines()
    #derniere ligne
    ligne = fichier[-1]
    tLigne = ligne.split()

    Tbulk=tLigne[1]

    ficTbulk.close()

    return Tbulk


def Re(mu,rho,D,U):
    nu = float(mu)/float(rho)
    reh = float(U)*float(D)/nu
    return reh

def Utheo(U,reh):
    lamb = 0.316*math.pow(reh,-0.25)
    utau = float(U)*math.sqrt(float(lamb)/8.00)
    return utau

def Uerr(Utheo,Utrio):
    Uerror = float(abs(Utheo-Utrio)/Utheo)*100
    return Uerror

def DHtheo(Utheo,D):
    DHth = (4*math.pow(float(Utheo),2))/float(D)
    return DHth

def DHtrio(Utrio,D):
    DHtrio=math.pow(Utrio,2)*4/float(D)
    return DHtrio

def DHtrio2():
    # ouverture des fichiers
    nomFic = os.popen('ls -rt *Pressure_Gradient_pb_sortie').readlines()[-1]
    nomFic = nomFic[:len(nomFic)-1] # Suppress /n

    fic = open(nomFic,'r')

    # lecture de ligne -> entetes
    fichier = fic.readlines()

    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    DHtrio=float(tLigne[1])
    fic.close()

    return DHtrio


def DHerr(DeltaHTheo,DeltaHTrio):
    DeltaHError = float(abs(DeltaHTheo-DeltaHTrio)/DeltaHTheo)*100
    return DeltaHError

def ecritureFichier(D,mu,lamb,Cp,rho,U,Ucalc,reh,rehCalc,yPlus,uPlus,Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2,DRey):

    nu = float(mu)/float(rho)

    #ecriture du fichier pour les vitesses et reynolds
    nomFicRe = 'ligneReynolds.dat'
    fichier = open(nomFicRe, 'w')
    fichier.write('%18.8f %18.8f %18.1f %18.1f %18.2f\n' % (U,Ucalc,reh,rehCalc,DRey))
    fichier.close()

    #ecriture du fichier pour la ligne du tableau
    nomFicLigne = 'ligneTableau.dat'
    fichier = open(nomFicLigne, 'w')
    fichier.write('%18.8f %18.8f %18.2f %18.8f %18.8f %18.2f %18.8f %18.2f\n' % (Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2))
    fichier.close()

    #ecriture du fichier pour le premier point
    nomFicFirstPoint = 'first_Point.dat'
    fichier = open(nomFicFirstPoint, 'w')
    fichier.write('%18.8f %18.8f\n' % (yPlus,uPlus))
    fichier.close()

    #ecriture du fichier pour la courbe Axial velocity
    nomFicRead = 'Conduite_SONDE_VP.coupe'
    nomFic = 'axial_velocity.dat'
    ficRead = open(nomFicRead, 'r')
    #1ere ligne vide
    ligne=ficRead.readline()
    fichier = open(nomFic, 'w')
    fin = False
    while not fin:
        ligne = ficRead.readline()
        if not ligne:
            fin=True
        else:
            tLigne = ligne.split()
            y = tLigne[0]
            u = tLigne[-1]
            yAdim = (float(y)+1)*float(Utr)/nu
            uAdim = float(u)/float(Utr)
            fichier.write('%18.8f %18.8f\n' % (yAdim,uAdim))
    fichier.close()
    ficRead.close()

def temperatureProfile(Tw,Qw,mu,rho,Cp,Utr):

    nu = float(mu)/float(rho)
    nomFicRead = 'Conduite_SONDE_TPN.coupe'
    nomFic = 'temperature_profile.dat'
    ficRead = open(nomFicRead,'r')
    fichier = open(nomFic,'w')
    #1ere ligne vide
    ligne=ficRead.readline()
    fin = False
    while not fin:
        ligne = ficRead.readline()
        if not ligne:
            fin=True
        else:
            tLigne = ligne.split()
            y = tLigne[0]
            t = tLigne[-1]
            yAdim = (float(y)+1)*float(Utr)/nu
            tAdim = (Tw-float(t))/(float(Qw)/(float(rho)*float(Cp)*float(Utr)))
            fichier.write('%18.8f %18.8f\n' % (yAdim,tAdim))
    fichier.close()
    ficRead.close()

def kader(mu,Cp,lamb,rho,Pr,Utr):

    nu = float(mu)/float(rho)

    nomFic = 'kader.dat'
    fichier = open(nomFic,'w')

    y=10
    while y<=2000:

        y=y+10
        beta = math.pow(3.85*math.pow(Pr,1.0/3.0)-1.3,2)+2.12*math.log(Pr)
        gamma = 0.01*math.pow(Pr*y,4)/(1+5*y*math.pow(Pr,3))
        tAdim = Pr*y*math.exp(-1*gamma)+(2.12*math.log(1+y)+beta)*math.exp(-1.0/gamma)
        fichier.write('%18.8f %18.8f\n' % (y,tAdim))

    fichier.close()


def metalLiquide(Pr,Utr):

    nu = float(mu)/float(rho)

    nomFic = 'metal_liquide.dat'
    fichier = open(nomFic,'w')

    y=10
    while y<=2000:

        y=y+10
        tAdim = (1.0/0.3)*math.log(1+0.3*Pr*y)
        fichier.write('%18.8f %18.8f\n' % (y,tAdim))

    fichier.close()


def PryPlus(Pr,Utr):
    nu = float(mu)/float(rho)

    nomFic = 'pryPlus.dat'
    fichier = open(nomFic,'w')

    y=10
    while y<=2000:

        y=y+10
        tAdim = Pr*y
        fichier.write('%18.8f %18.8f\n' % (y,tAdim))

    fichier.close()

def nusselt(TwTotal,Qw,D,lamb,Tbulk,reh,Pr):
    nomFic = 'nusselt.dat'
    fichier = open(nomFic,'w')
    Nu=float(Qw)*float(D)/float(lamb)/(TwTotal-float(Tbulk))

    fichier.write('%18.8f %18.8f %18.8f\n' % (Nu,reh,Pr))
    fichier.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    #recuperation des donnees
    D,mu,lamb,Cp,rho,U,Qw,Z = properties()

    Pr = float(mu)*float(Cp)/float(lamb)

    reh = Re(mu,rho,D,U)
    Ucalc = UbulkCalc(rho,D)
    rehCalc = Re(mu,rho,D,Ucalc)
    y = yPlus()
    u = uPlus()
    Uth = Utheo(U,reh)
    Uthcalc = Utheo(Ucalc,rehCalc)
    Utr = Utrio() #u*
    Uerror = Uerr(Uth,Utr)
    DeltaHTheo = DHtheo(Uth,D)
    DeltaHTrio = DHtrio(Utr,D)
    DeltaHTrio2 = DHtrio2()
    DeltaHError = DHerr(DeltaHTheo,DeltaHTrio)
    DeltaHError2 = DHerr(DeltaHTheo,DeltaHTrio2)
    DRey = DHerr(reh,rehCalc)

    Tw = Tparoi(Z)

    TwTotal = Ttotal()
    Tbulk = Tbulk()

    #ecriture des fichiers gnuplot et tableaux
    ecritureFichier(D,mu,lamb,Cp,rho,U,Ucalc,reh,rehCalc,y,u,Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2,DRey)

    #ecriture du fichier gnuplot pour T+=f(y+)
    temperatureProfile(Tw,Qw,mu,rho,Cp,Utr)

    #ecriture du fichier gnuplot pour Kader
    kader(mu,Cp,lamb,rho,Pr,Utr)

    #ecriture du fichier gnuplot pour la loi metal liquide
    metalLiquide(Pr,Utr)

    #ecriture du fichier gnuplot pour la courbe Pry+
    PryPlus(Pr,Utr)

    #ecriture du nusselt
    nusselt(TwTotal,Qw,D,lamb,Tbulk,reh,Pr)
