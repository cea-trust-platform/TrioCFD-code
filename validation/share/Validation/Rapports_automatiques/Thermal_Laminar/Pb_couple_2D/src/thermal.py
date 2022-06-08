import os,os.path, sys, math

def getflux(ADC, APT, FCC, FDC):
#
    fadc = open(ADC, 'r')
    fapt = open(APT, 'r')
    ffcc = open(FCC, 'r')
    ffdc = open(FDC, 'r')

    ladc = fadc.readlines()
    lapt = fapt.readlines()
    lfcc = ffcc.readlines()
    lfdc = ffdc.readlines()
    fadc.close()
    fapt.close()
    ffcc.close()
    ffdc.close()
#
    fic = open('thermal.dat', 'w')
#adc
    liadc = ladc[-1].strip()
    tLiadc = liadc.split()
    adc6 = float(tLiadc[5])
#apt
    liapt = lapt[-1].strip()
    tLiapt = liapt.split()
    apt2 = float(tLiapt[1])
#fcc
    lifcc = lfcc[-1].strip()
    tLifcc = lifcc.split()
    fcc6 = float(tLifcc[5])
#fdc
    lifdc = lfdc[-1].strip()
    tLifdc = lifdc.split()
    fdc2 = float(tLifdc[1])
    fdc3 = float(tLifdc[2])
#
    flux = 180 + (fcc6+fdc2+fdc3)
    fic.write('%18.1f %18.1f %18.1f %18.1f %18.1f %18.f\n' % (adc6,apt2,fcc6,fdc2,fdc3,flux))
#
if __name__=='__main__':
    #recuperation des parametres passes en ligne de commande
    args = sys.argv
    if len(args)!=5:
        print('Erreur sur le nb d\'arguments fournis : Usage\npython calculerFrequence nomFichiers.out')
        sys.exit()

    ADC = args[1]
    APT = args[2]
    FCC = args[3]
    FDC = args[4]
    getflux(ADC, APT, FCC, FDC)
