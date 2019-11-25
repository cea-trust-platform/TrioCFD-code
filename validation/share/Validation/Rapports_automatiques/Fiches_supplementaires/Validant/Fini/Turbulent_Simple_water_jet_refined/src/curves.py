import os,os.path, sys, math
###
###
### Fonction de recuperation des parametres issus di fichier data :
# (fichier genere par propertiesGeometry.py)
#  -vitesse imposee
#  -concentration imposee
#  -diametre du jet
# ...
# cf propertiesGeometry.py
def properties(propertiesFile):
    # ouverture du fichier de proprietes
    fic = open(propertiesFile,'r')
    lignes = fic.readlines()
    fic.close()
    l = lignes [0]
    tLi = l.split()
    vitesseImposee = float(tLi[3])
    concentrationImposee = float(tLi[4])
    diametreJet = float(tLi[5])
    velocityDecayCte = float(tLi[6])
    velocityOrigin = float(tLi[7])
    concentrationDecayCte = float(tLi[8])
    concentrationOrigin = float(tLi[9])
    return vitesseImposee,concentrationImposee, diametreJet, velocityDecayCte,velocityOrigin, concentrationDecayCte,concentrationOrigin
##
#
# dimensionless concentration
# Methode de generation du fichier cc.dat
#Ce fichier contient les valeurs de la sonde de concentration, pour z dans [0.003,0.15] :
# zAdim concentrationAdim z concentration
def adimensionConcentration(fc,concentrationImposee,diametreJet):
    p = 0
    ffc = open('cc.dat','w')
    fcc = open(fc, 'r')
    lignes = fcc.readlines()
    for ligne in lignes:
        p +=1
        if  p>1:
            ligne = ligne.strip()
            tLigne = ligne.split()
            zz = float(tLigne[0])
            if (p==2):
                c_centre=float(tLigne[1])
                #print 'c_centre=', c_centre
            if ((zz>=0.003)&(zz<=0.15)):
                c1 = zz/diametreJet
                c12 = zz
                c2 = float(tLigne[1])/c_centre
                #c2 = float(tLigne[1])/concentrationImposee
                c21 = float(tLigne[1])
                ffc.write('%18.3f %18.3f %18.3f %18.3f\n' % (c1,c2,c12,c21))
    fcc.close()
#
# dimensionless velocity
#Ce fichier contient les valeurs de la sonde de vitesse, pour z dans [0.003,0.15] :
# zAdim vitesseAdim z vitesse
def adimensionVelocity(sondeVitesseFile,vitesseImposee,diametreJet):
    p = 0
    ffuc = open('uc.dat','w')
    fuc = open(sondeVitesseFile, 'r')
    lignes = fuc.readlines()
    fuc.close()
    for ligne in lignes:
        p +=1
        if  p>1:
            ligne = ligne.strip()
            tLigne = ligne.split()
            zz = float(tLigne[0])
            if (p==2):
                v_centre=float(tLigne[1])
                #print 'v_centre=', v_centre
            if ((zz>=0.003)&(zz<=0.15)):
                c1 = zz/diametreJet
                c12 = zz
                c2 = float(tLigne[1])/v_centre
                #c2 = float(tLigne[1])/vitesseImposee
                #print 'vitesseImposee=%f' % (vitesseImposee)
                c21 = float(tLigne[1])
                ffuc.write('%18.4f %18.4f %18.4f %18.4f\n' % (c1,c2,c12,c21))
#
# Momentum jet width
# Calcul de la largeur du jet, en fonction de la vitesse ou de la concentration
# La vitesse centrale theorique Uc(z) est donnee par (Uj etant la vitesse d'injection, d le diametre d'injection, C1 la constante de reduction de la vitesse) :
#    Uc(z) = Uj . C1/( (z-z0)/d )
# la largeur b du jet correspond au rayon pour lequel la vitesse longitudinale devient inferieure a Uc(z)/e
# en pratique, dans ce calcul, on prend Uc calculee
# (description ci-dessus valable en remplacant "vitesse" par "concentration"
#Methode de calcul de b, a z donne
# type = U (calcul en vitesse) ou C (calcul en concentration)
def calculerBz(z, type, valJet, valDecayConstante, valOrigin, diametreJet):
    b = -1
    valadim = []
    if type in ['U', 'C']:
        #Recupere le nom de la sonde
        nomFic = 'test_SONDE_%sZ%5.3f.coupe' % (type, z)
        if os.path.exists(nomFic):
            #recupere les donnees du fichier
            f = open(nomFic, 'r')
            lignes = f.readlines()
            f.close()
            #calcule la grandeur de reference Uc ou Cc, et cette valeur/e
            #print 'type=',type, '  z=',z, ' valjet=',valJet, ' C=', valDecayConstante, ' valorigin=',valOrigin, ' diam=',diametreJet
            valC = valJet * valDecayConstante / ( ( z - valOrigin*diametreJet ) / diametreJet) # / math.exp(1)
            valC_e = valC  / math.exp(1)
            #print ' => Uc=',valC, '   Uc/e=',valC_e
            rprec = 0.
            valprec = 0.
            b = 0.
            valCCalc = -1000000000.
            valC_eCalc = -1000000000.
            #balaye le rayon
            for ligne in lignes:
                ligne = ligne.strip()
                if len(ligne)>0:
                    ligneT = ligne.split()
                    r = float(ligneT[0])
                    val = float(ligneT[1])
                    #memorise la 1e valeur de la vitesse = Uc calculee
                    if valCCalc<-99999.:
                        valCCalc = val
                        valC_eCalc = valCCalc / math.exp(1)
                    if val<=valC_eCalc:
                        #on a traverse le point tq U/Uc=1/e  => recupere b par interpolation lineaire
                        b = r + ( (valC_eCalc - val) / (valprec - val) ) * (rprec - r)
                        #print '       => b=',b
                        break
                    rprec = r
                    valprec = val
            #print '  valCTh=',valC, '   valCCalc=',valCCalc
            for ligne in lignes:
                ligne = ligne.strip()
                if len(ligne)>0:
                    ligneT = ligne.split()
                    r = float(ligneT[0])
                    val = float(ligneT[1])
                    #calcul des valeurs adimensionnees
                    valadim.append([r/b, val/valCCalc])
                    #print '    r=',r, '  val=', val, '  ==>> adim=',r/b,val/valCCalc
        else:
            print 'Erreur : fichier %s non trouve !' % (nomFic)
    else:
        print 'Erreur : type %s non reconnu (U ou C)' % (type)

    return b, valadim


def calculerJetWidth(nomFic, type, valImposee, valDecayCte, valOrigin, diametreJet):
    zz = [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15]
    fb = open(nomFic, 'w')
    for z in zz[1:]:
        b, valadim = calculerBz(z, type, valImposee, valDecayCte, valOrigin, diametreJet)
        #ecriture de b
        fb.write('%18.5f %18.5f %18.5f\n' % (z/diametreJet, b/diametreJet, b))
        #ecriture des grandeurs adim dans un fichier
        nomFicRadialVal = 'radial%s_%5.3f.dat' % (type, z)
        fr = open(nomFicRadialVal, 'w')
        for el in valadim:
            fr.write('%18.3f %18.3f\n' % (el[0], el[1]))
        fr.close()
    fb.close()

def calculerMomentumJetWidth(vitesseImposee, velocityDecayCte,velocityOrigin, diametreJet):
    nomFicRes = 'bvalues.dat'
    calculerJetWidth(nomFicRes, 'U', vitesseImposee, velocityDecayCte,velocityOrigin, diametreJet)
def calculerConcentrationJetWidth(concentrationImposee, concentrationDecayCte,concentrationOrigin, diametreJet):
    nomFicRes = 'btvalues.dat'
    calculerJetWidth(nomFicRes, 'C', concentrationImposee, concentrationDecayCte,concentrationOrigin, diametreJet)


if __name__=='__main__':
    #recuperation des parametres passes en ligne de commande
    args = sys.argv
    if len(args)!=4:
        print 'Erreur sur le nb d\'arguments fournis : Usage\npython sondeVitesse.coupe sondeConcentration.coupe propertiesGeometry.dat'
        sys.exit()
#
    sondeVitesseFile = args[1]
    fc = args[2]
    propertiesFile = args[3]
    #recupere les vitesse & concentration imposees
    #et le diametre du jet dans le fichier propertiesGeometry.dat
    vitesseImposee, concentrationImposee, diametreJet, velocityDecayCte,velocityOrigin, concentrationDecayCte,concentrationOrigin = properties(propertiesFile)
    #adimensionnalise les vitesse & concentration
    adimensionConcentration(fc,concentrationImposee,diametreJet)
    adimensionVelocity(sondeVitesseFile,vitesseImposee,diametreJet)

    #Calcule les largeurs du jet
    #Momentum jet width
    calculerMomentumJetWidth(vitesseImposee, velocityDecayCte,velocityOrigin, diametreJet)
    #Concentration jet width
    calculerConcentrationJetWidth(concentrationImposee, concentrationDecayCte,concentrationOrigin, diametreJet)
