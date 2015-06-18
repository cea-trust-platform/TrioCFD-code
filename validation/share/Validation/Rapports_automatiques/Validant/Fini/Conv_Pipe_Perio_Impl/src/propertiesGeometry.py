#Script de recuperation des propretes physiques et de la geometrie

import optparse
import math
import os

def properties(nomFic):
	# ouverture des fichiers
	fic = open(nomFic,'r')
	
	chaines = ["mu Champ_Uniforme",
		   "lambda Champ_Uniforme",
		   "beta_th Champ_Uniforme",
		   "Cp Champ_Uniforme",
		   "rho Champ_Uniforme",
		   "vitesse champ_fonc_xyz",
		   "paroi   paroi_flux_impose Champ_Front_Uniforme",
		   "sonde_tpn "] # Texte a rechercher
	
	for ligne in fic:
	    for chaine in chaines:
	        if chaine in ligne:
	            tLigne = ligne.split()
		    if chaine=="mu Champ_Uniforme":
			mu=float(tLigne[-1])
		    if chaine=="lambda Champ_Uniforme":
			lamb=float(tLigne[-1])
		    if chaine=="Cp Champ_Uniforme":
			Cp=float(tLigne[-1])
		    if chaine=="rho Champ_Uniforme":
		    	rho=float(tLigne[-1])
		    if chaine=="vitesse champ_fonc_xyz":
		   	formule = tLigne[-1]
			tFormule = formule.split('*')
			U=float(tFormule[0])/2.000
		    if chaine=="paroi   paroi_flux_impose Champ_Front_Uniforme":
		    	Qw=float(tLigne[-1])
		    if chaine=="sonde_tpn ":
		    	Z=float(tLigne[-1])
		    	D=2.*float(tLigne[9])
	fic.close()
	return D,mu,lamb,Cp,rho,U,Qw,Z


def ecritureFichier(D,mu,lamb,Cp,rho,U,Qw,Z):
	#ecriture du fichier
	nomFic = 'propertiesGeometry.dat'
	fichier = open(nomFic, 'w')
	fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.8f %18.2f\n' % (D,mu,lamb,Cp,rho,U,Qw,Z))
	fichier.close()

if __name__ == '__main__':
	
	#recuperation du fichier data
	#derniere ligne du ls
	ficLS = os.popen('ls *.data')
	lignes = ficLS.readlines()
	Ligne = lignes[0]
	#suppression du \n en fin de nom
	nomFic = Ligne[:len(Ligne)-1]
	D,mu,lamb,Cp,rho,U,Qw,Z = properties(nomFic)
	
	#ecriture du fichier
	ecritureFichier(D,mu,lamb,Cp,rho,U,Qw,Z)

