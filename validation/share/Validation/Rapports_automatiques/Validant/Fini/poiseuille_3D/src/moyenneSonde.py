#Script de calcul de la moyenne des profils de vitesse en horizontal et vertical

import optparse


def moyenneSonde(ficH,ficV):
	# ouverture des fichiers
	ficHOR = open(ficH,'r')
	ficVERT = open(ficV,'r')
	# lecture de ligne -> entetes
	ligneH = ficHOR.readline()
	ligneV = ficVERT.readline()

	#recuperation des valeurs
	Temps=[]
	VZ=[]
	fin=False
	nbPoints = 0
	while not fin:
		ligneH = ficHOR.readline()
		ligneV = ficVERT.readline()
		if not ligneH or not ligneV:
			fin=True
		tLigneH = ligneH.split()
		tLigneV = ligneV.split()
		if len(tLigneH)>0 and len(tLigneV)>0:
			tps=float(tLigneH[0])
			Temps.append(tps)
			zh=float(tLigneH[3])
			zv=float(tLigneH[3])
			VZ.append((zh+zv)/2)
			nbPoints += 1
	ficHOR.close()
	ficVERT.close()
	return Temps, VZ, nbPoints

def ecritureFichier(X,Y,nbpoints):
	#ecriture du fichier pour gnuplot
	fichier = open('./moyenneSonde.plot', 'w')
	i = 0
	while i<nbpoints:
		fichier.write('%18.8f %18.8f\n' % (X[i],Y[i]))
		i += 1
	fichier.close()


if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	#recuperation des moyennes des valeurs et des pas de temps
	Temps,VZ,nbpoints = moyenneSonde(args[0],args[1])
	
	#ecriture du fichier gnuplot
	ecritureFichier(Temps,VZ,nbpoints)

