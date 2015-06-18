#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties():
	# ouverture des fichiers
	nomFic = '../Calcul.data'
			
	fic = open(nomFic,'r')
	
	chaines = ["mu Champ_Uniforme",
	           "rho Champ_Uniforme"] # Texte a rechercher

	for ligne in fic:
	    for chaine in chaines:
	        if chaine in ligne:
	            tLigne = ligne.split()
		    
		    if chaine=="mu Champ_Uniforme":
			mu=float(tLigne[3])
		    if chaine=="rho Champ_Uniforme":
			rho=float(tLigne[3])

	fic.close()
	return mu,rho


def ecritureFichier(mu,rho):
	#ecriture du fichier
	nomFic = 'properties.dat'
	fichier = open(nomFic, 'w')
	fichier.write('%18.8f %18.8f\n' % (mu,rho))
	fichier.close()

if __name__ == '__main__':
	
	mu,rho = properties()
	
	#ecriture du fichier
	ecritureFichier(mu,rho)

