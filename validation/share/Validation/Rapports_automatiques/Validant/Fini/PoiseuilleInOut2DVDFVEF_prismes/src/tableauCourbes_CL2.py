#Script de construction des tableaux 
#Pour cette fiche, fonctionne pour le cas ou Ventree et Psortie sont imposes a 0

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
	U0=float(tLigne[4])
	fic.close()
	return mu,rho,h,L,U0

#CALCUL DU REYNOLDS GRACE A LA VITESSE MOYENNE DANS LE CANAL. FORMULE ANALYTIQUE
def Re(U0,mu,h,rho):    
	
	reh=4*rho*U0*h/mu	
	return reh

def Utheo(reh,U0):

	utau = float(U0)*math.sqrt(12/reh)
	return utau

def Utrio():
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
	utau=float(tLigne[-1])
	ficUstar.close()
	
	return utau
	
def Uerr(Utheo,Utrio):
	Uerror = float(abs(Utheo-Utrio)/Utheo)*100
	return Uerror

def Ptheo(Utheo,L,h,rho):
	Pth = float(Utheo*Utheo*rho*L/h)
	return Pth

def Ptrio(L):
	# ouverture des fichiers
	nomFic = 'test_SEG_P.son'
			
	fic = open(nomFic,'r')
	
	# lecture de ligne -> entetes
	fichier = fic.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	ptrio=float(tLigne[1])-float(tLigne[-1])
	fic.close()
	
	return ptrio
	
def Perr(Ptheo,Ptrio):
	Perror = float(abs(Ptheo-Ptrio)/Ptheo)*100
	return Perror

def ecritureFichier(U,Re,Utheo,Utrio,UError,Ptheo,Ptrio,PError):
	#ecriture du fichier pour la ligne du tableau
	nomFicLigne = 'ligneTableau.dat'
	fichier = open(nomFicLigne, 'w')
	fichier.write('%18.4f %18.4f %18.4f %18.4f %18.4f %18.4f %18.4f %18.4f\n' % (U,Re,Utheo,Utrio,UError,Ptheo,Ptrio,PError))
	fichier.close()

def copieFichierMoyennesSpatiales(nomFic):
	#ecriture du fichier
	ficRead = open(nomFic,'r')
	ficWrite = open('Velocity_profile.dat', 'w')

	fin=False
	i=0
	while not fin:
		ligne = ficRead.readline()
		if not ligne:
			fin=True
		else:
			if ligne[0]!="#" and ligne[0]!=" " and ligne[0]!="":
				ficWrite.write('%s' % (ligne))
	
	ficWrite.close()
	ficRead.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	mu,rho,h,L,U0 = properties()
	
	reh = Re(U0,mu,h,rho)
	Uth = Utheo(reh,U0)
	Utr = Utrio()
	Uerror = Uerr(Uth,Utr)
	Pth = Ptheo(Uth,L,h,rho)
	Ptr = Ptrio(L)
	Perror = Perr(Pth,Ptr)
	
	#ecriture des fichiers gnuplot et tableaux
	ecritureFichier(U0,reh,Uth,Utr,Uerror,Pth,Ptr,Perror)
	
	#derniere ligne du ls
	ficLS = os.popen('ls -rt Moyennes_spatiales_vitesse_rho_mu_*')
	lignes = ficLS.readlines()
	derLigne = lignes[-1]
	#suppression du \n en fin de nom
	nomFic = derLigne[:len(derLigne)-1]
	
	#copie du fichier MoyennesSpatiales* dans un fichier generique VelocityProfile
	copieFichierMoyennesSpatiales(nomFic)
