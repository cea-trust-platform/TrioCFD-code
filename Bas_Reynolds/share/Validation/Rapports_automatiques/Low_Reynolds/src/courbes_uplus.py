#Script de recuperation des donnees u* et calcul de U+ et y+

import optparse
import math
import os

def properties():
	# ouverture des fichiers
	nomFic = 'propertiesGeometry.dat'
			
	fic = open(nomFic,'r')
	
	# lecture de ligne -> entetes
	fichier = fic.readlines()
	
	ligne = fichier [0]
	tLigne = ligne.split()
	mu=float(tLigne[0])
	rho=float(tLigne[1])
	h=float(tLigne[2])
	L=float(tLigne[3])
	
	fic.close()
	return mu,rho,h,L

#fonction de recuperation de la colonne u+
def utau():
	# ouverture des fichiers
	nomFicUstar = 'u_tau.dat'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	utau=float(tLigne[3])
	ficUstar.close()
	print "utau",utau
	return utau


def ecritureFichier(utau,nomFicData,rho,mu):
	#ecriture du fichier pour gnuplot
	nomFic = 'courbe_uplus.dat'
	ficRead = open(nomFicData,'r')
	fichier = open(nomFic, 'w')
	
	# lecture de ligne -> entetes
	i=0
	while i<17:
		ligne = ficRead.readline()
		i=i+1

	#recuperation des valeurs
	fin=False
	while not fin:
		ligne = ficRead.readline()
		if not ligne:
			fin=True
		tLigne = ligne.split()
		if len(tLigne)>0:
			y=float(tLigne[0])
			yPlus=y*rho*float(utau)/mu
			u=float(tLigne[1])
			uPlus=u/utau
			fichier.write(' %18.8f   %18.8f   %18.8f  %18.8f\n' % (y,yPlus,u,uPlus))
	
	fichier.close()
	ficRead.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	mu,rho,h,L = properties()
	
	#recuperation des donnees
	utau = utau()
	
	#recuperation du fichier avec le dernier temps dans le nom
	#derniere ligne du ls
	ficLS = os.popen('ls -rt '+args[0]+'*')
	lignes = ficLS.readlines()
	derLigne = lignes[-1]
	#suppression du \n en fin de nom
	nomFic = derLigne[:len(derLigne)-1]
	
	#ecriture du fichier gnuplot
	ecritureFichier(utau,nomFic,rho,mu)

