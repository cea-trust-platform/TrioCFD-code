import sys, os.path

#fonction de recuperation de la valeur u_tau
def getUTau(dir):
    # ouverture des fichiers
    nomFicUstar = os.path.join(dir, 'test_pb_Ustar.face')

    ficUstar = open(nomFicUstar,'r')

    # lecture de ligne -> entetes
    fichier = ficUstar.readlines()

    #en periodique on prend la valeur moyenne calculee en fin de fichier
    #tant que la derniere ligne est vide
    while fichier[-1]=="" or fichier[-1]=="\n":
        del fichier [-1]

    ligne = fichier [-1]
    tLigne = ligne.split()
    utau=float(tLigne[5])
    ficUstar.close()

    return utau

if __name__=='__main__':
    #recuperation des parametres passes en ligne de commande : les repertoires ou recuperer les u_tau
    #dans les fichiers test_pb_Ustar.face
    args = sys.argv
    res = ''
    utauTH = float(args[1])
    for dir in args[2:]:
        if os.path.isdir(dir):
            utau = getUTau(dir)
            res += '  %10.2f %10.2f   ' % (utau, abs(utau-utauTH)/utauTH*100.)
    res += '   %10.2f' % utauTH
    fic = open('utau.dat', 'w')
    fic.write('%s\n' % res)
    fic.close()
