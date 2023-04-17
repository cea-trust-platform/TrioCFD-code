import sys, os.path

#fonction de recuperation de la colonne y+
def getYPlus(dir):
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
    y=float(tLigne[3])
    ficUstar.close()

    return y

if __name__=='__main__':
    #recuperation des parametres passes en ligne de commande : les repertoires ou recuperer les y+
    #dans les fichiers test_pb_Ustar.face
    args = sys.argv
    res = ''
    for dir in args[1:]:
        if os.path.isdir(dir):
            y = getYPlus(dir)
            res += '  %18.2f' % y
    fic = open('yplus.dat', 'w')
    fic.write('%s\n' % res)
    fic.close()
