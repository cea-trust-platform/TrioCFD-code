#Script de construction des tableaux et des courbes

import optparse
import math

def ecritureFichier():

    nomFicEFstab = 'EF_stab/nusselt.dat'
    nomFicPrtStd = 'Prt_std/nusselt.dat'
    nomFicMuscl = 'Muscl/nusselt.dat'

    nomFic = 'tableau_nusselt.dat'

    ficEFstab = open(nomFicEFstab,'r')
    ficPrtStd = open(nomFicPrtStd,'r')
    ficMuscl = open(nomFicMuscl,'r')

    fichier = open(nomFic,'w')

    ligneEFstab=ficEFstab.readline()
    tLigneEFstab=ligneEFstab.split()
    NuEFstab=tLigneEFstab[0]
    reh=tLigneEFstab[1]
    Pr=tLigneEFstab[2]

    lignePrtStd=ficPrtStd.readline()
    tLignePrtStd=lignePrtStd.split()
    NuPrtStd=tLignePrtStd[0]

    ligneMuscl=ficMuscl.readline()
    tLigneMuscl=ligneMuscl.split()
    NuMuscl=tLigneMuscl[0]

    NuH=6.3+0.0167*math.pow(float(reh),0.85)*math.pow(float(Pr),0.93)
    varEFstab=float(abs(NuH-float(NuEFstab))/NuH)*100
    varPrtStd=float(abs(NuH-float(NuPrtStd))/NuH)*100
    varMuscl=float(abs(NuH-float(NuMuscl))/NuH)*100

    fichier.write('%18.2f %18.2f %18.2f %18.2f %18.2f %18.2f %18.2f\n' % (NuH,float(NuEFstab),varEFstab,float(NuPrtStd),varPrtStd,float(NuMuscl),varMuscl))

    fichier.close()

    ficEFstab.close()
    ficPrtStd.close()
    ficMuscl.close()


if __name__ == '__main__':

    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    #ecriture du fichier tableau Nusselt
    ecritureFichier()
