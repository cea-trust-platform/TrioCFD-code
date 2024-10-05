import sys, os.path, math
import propertiesGeometry


def getValeurMax(nomFichier, colonne):
    #renvoie la valeur max de la colonne 'colonne', avec la coordonnee correspondante (1e colonne)
    if os.path.isfile(nomFichier):
        #recupere les donnees du fichier
        f = open(nomFichier, 'r')
        lignes = f.readlines()
        f.close()
        #balaye les donnees pour trouver al valeur max
        xmax = -1
        valmax = -1
        for ind, ligne in enumerate(lignes):
            ligne = ligne.strip()
            if len(ligne)>0:
                tabLigne = ligne.split()
                try:
                    val = float(tabLigne[colonne])
                except IndexError:
                    print('Error : line (%d : %s) has no %d values !' % (ind, ligne, colonne))
                    sys.exit()
                except ValueError:
                    print('Error : line (%d : %s) does not contain float at index  %d !' % (ind, ligne, colonne))
                    sys.exit()
                if val>valmax:
                    valmax = val
                    xmax = float(tabLigne[0])
        return xmax, valmax
    else:
        print('Error: file %s not found!' % (nomFichier))
        sys.exit()

if __name__ == '__main__':
    #recuperation des parametres passes en ligne de commande
    args = sys.argv
    if len(args)!=4:
        print(r'''Error: wrong number of passed arguments: Usage
python extraitMax.py LongueurCarac Tparoi Tfluide''')
        sys.exit()

    L = float(args[1])
    Tparoi = float(args[2])
    Tfluide = float(args[3])

    #couples des cas [nom de fichier, colonne ou chercher le max]
    cas = [['test_SONDE_VIT1.coupe', 1],
            ['test_SONDE_VIT0.coupe', 2]]
    ref = [
            ['Vahl_Davis', 64.63, 216.36],
            ['Gresho', 64.593, 220.64],
            ['Winter', 63.9, 222.],
            ]

    properties = propertiesGeometry.getPropertiesFromdat()
    print('Physical properties= %s' % (properties))

    f = open('vitesseMaxAdim.dat', 'w')
    #f.write('#xmax VitAdim Rayleigh Prandtl refVahlDavis errVahlsDavis refGresho errGresho refWinter errWinter\n')
    for nomFichier, colonne in cas:
        xmax, valmax = getValeurMax(nomFichier, colonne)
        print('Max value %f not found at %f' % (valmax, xmax))

        #adimensionnalisation
        xadim = xmax / L
        Prandtl = propertiesGeometry.nombreDePrandtl(properties['mu'], properties['lambda'], properties['cp'])
        print('Prandtl=%f' % (Prandtl))
        Rayleigh = propertiesGeometry.nombreDeRayleigh(properties['g'], properties['beta_th'], properties['rho'], properties['mu'], properties['alpha'], Tparoi, Tfluide, L)
        print('Rayleigh=%f' % (Rayleigh))
        valadim = valmax * L * math.sqrt(Rayleigh * Prandtl)
        f.write('%18.8f %18.8f %18.8f %18.8f' % (valadim, xmax, Rayleigh, Prandtl))
        for valsref in ref:
            valref = valsref[colonne]
            maxx = max(valadim, valref)
            err = 100 * abs(valadim - valref) / maxx
            print(r'''Probe %s Max dimensionless value = %f found at %f
Relative error (%s : %f) %f''' % (nomFichier, valadim, xadim, valsref[0], valref, err))
            f.write(' %18.8f %18.8f' % (valref, err))
        f.write('\n')
    f.close()
