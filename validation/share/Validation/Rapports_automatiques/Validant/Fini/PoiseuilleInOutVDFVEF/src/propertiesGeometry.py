#Script de recuperation des proprietes physiques et de la geometrie

import os
import optparse
import sys
import math

def readPropertiesData():
	# ouverture des fichiers
	nomFic = 'test.data'
			
	fic = open(nomFic,'r')
	
	chaines = ["Longueurs",
        	   "mu Champ_Uniforme",
	           "rho Champ_Uniforme",
		   "Entree frontiere_ouverte_pression_imposee Champ_Front_Uniforme",
	           "Sortie frontiere_ouverte_pression_imposee Champ_Front_Uniforme"] # Texte a rechercher

	for ligne in fic:
	    for chaine in chaines:
	        if chaine in ligne:
	            tLigne = ligne.split()
		    if chaine=="Longueurs":
		    	h=float(tLigne[2])
			L=float(tLigne[1])
		    if chaine=="mu Champ_Uniforme":
			mu=float(tLigne[3])
		    if chaine=="rho Champ_Uniforme":
			rho=float(tLigne[3])
		    if chaine=="Entree frontiere_ouverte_pression_imposee Champ_Front_Uniforme":
			Pin=float(tLigne[4])
		    	print 'Pin:%f' %Pin 
         	    if chaine=="Sortie frontiere_ouverte_pression_imposee Champ_Front_Uniforme":
			Pout=float(tLigne[4])

	fic.close()
	return mu,rho,h,L,Pin,Pout
	return properties


def ecritureFichier(mu,rho,h,L,Pin,Pout):
	#ecriture du fichier
	nomFic = 'propertiesGeometry.dat'
	fichier = open(nomFic, 'w')
	fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.8f\n' % (mu,rho,h,L,Pin,Pout))
	fichier.close()




if __name__ == '__main__':
	
	
	mu,rho,h,L,Pin,Pout = readPropertiesData()
	
	#ecriture du fichier
	ecritureFichier(mu,rho,h,L,Pin,Pout)
	

