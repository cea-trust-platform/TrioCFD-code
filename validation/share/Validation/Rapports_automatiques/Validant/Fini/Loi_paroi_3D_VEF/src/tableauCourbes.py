#Script de construction des tableaux et des courbes

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
	Lx=float(tLigne[0])
	Ly=float(tLigne[1])
	Lz=float(tLigne[2])
	mu=float(tLigne[3])
	rho=float(tLigne[4])
	U=float(tLigne[5])	
	
	fic.close()
	return Lx,Ly,Lz,mu,rho,U

def Utrio():
	# ouverture des fichiers
	nomFicUstar = 'u_tau.dat'
			
	ficUstar = open(nomFicUstar,'r')
	
	#on prend seulement les lignes entre 20 et 35 secondes
	seconde=-1
	while seconde<20:
		ligne = ficUstar.readline()
		tLigne = ligne.split()
		tps = tLigne[0]
		try:
		 seconde =float(tps)
		except:
		 print "pas un temps",tps
		 pass

		
	#lecture de l ensemble des donnees restantes
	tFic = ficUstar.readlines()
	
	#tant que la derniere ligne est vide on la supprime
	while tFic[-1]=="" or tFic[-1]=="\n":
		del tFic [-1]
	
	#on ajoute la premiere ligne a 20 secondes
	tFic.append(ligne)

	i=len(tFic)
	cumul=0
	nbTemps=0
	
	while i>0:
		ligne=tFic[-i]
		tLigne=ligne.split()
		cumul=cumul+float(tLigne[1])
		nbTemps=nbTemps+1
		i=i-1

	utau=cumul/nbTemps
	ficUstar.close()
	
	return utau

def forcTrio():
	# ouverture des fichiers
	nomFic = os.popen('ls -rt *Pressure_Gradient_pb_periox').readlines()[-1]
	nomFic = nomFic[:len(nomFic)-1] # Suppress /n
			
	fic = open(nomFic,'r')
	
	#on prend seulement les lignes entre 20 et 35 secondes
	seconde=-1
	while seconde<20:
		ligne = fic.readline()
		tLigne = ligne.split()
		tps = tLigne[0]
		try:
		 seconde =float(tps)
		except:
		 print "pas un temps",tps
		 pass
		
	#lecture de l ensemble des donnees restantes
	tFic = fic.readlines()
	
	#tant que la derniere ligne est vide on la supprime
	while tFic[-1]=="" or tFic[-1]=="\n":
		del tFic [-1]
	
	#on ajoute la premiere ligne a 20 secondes
	tFic.append(ligne)

	i=len(tFic)
	cumul=0
	nbTemps=0
	
	while i>0:
		ligne=tFic[-i]
		tLigne=ligne.split()
		cumul=cumul+float(tLigne[1])
		nbTemps=nbTemps+1
		i=i-1

	forc=cumul/nbTemps
	fic.close()
	
	return forc

def Utheo(Re_tau,Ly,rho,mu):
	utau = float(Re_tau)*(float(Ly)/2.0)*float(mu)/float(rho)
	return utau

def forcTheo(Uth,Ly):
	forc = (2.0/float(Ly))*math.pow(Uth,2)
	return forc	

def Uerr(Utheo,Utrio):
	Uerror = float(abs(Utheo-Utrio)/Utheo)*100
	return Uerror

def ecritureFichier(Uth,forcTh,Utr,forcTr,Uerror,nu):
	
	#ecriture du fichier pour la ligne theorique du tableau
	nomFicLigne = 'ligneTheoTableau.dat'
	fichier = open(nomFicLigne, 'w')
	fichier.write('%18.8f %18.8f\n' % (Uth,forcTh))
	fichier.close()

	#ecriture du fichier pour la ligne trio_U du tableau
	nomFic = 'ligneTableau.dat'
	fichier = open(nomFic, 'w')
	fichier.write('%18.8f %18.2f %18.8f\n' % (Utr,Uerror,forcTr))
	fichier.close()
	 
	#ecriture du fichier pour la courbe adimensionnee	
	nomFicRead = 'Moyennes_spatiales_vitesse_rho_mu'
	nomFic = 'velocity_profile.dat'
	ficRead = open(nomFicRead, 'r')
	#lignes vides
	i=0
	while i<18:
		ligne=ficRead.readline()
		i=i+1
	
	
	fichier = open(nomFic, 'w')
	fin = False
	while not fin:
		ligne = ficRead.readline()
		if not ligne:
			fin=True
		else:
			if ligne!="" and ligne!="\n":
				tLigne = ligne.split()
				y = tLigne[0]
				u = tLigne[1]
				yAdim = (float(y))*float(Utr)/nu #(float(y)+1)*float(Utr)/nu
				uAdim = float(u)/float(Utr)
				fichier.write('%18.8f %18.8f\n' % (yAdim,uAdim))
	fichier.close()
	ficRead.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	#recuperation des donnees
	Lx,Ly,Lz,mu,rho,U = properties()
	
	nu = float(mu)/float(rho)
	
	Re_tau = float(rho)*float(U)*2*float(Ly)/float(mu)
	
	Uth = Utheo(Re_tau,Ly,rho,mu)
	
	forcTh = forcTheo(Uth,Ly)
	
	Utr = Utrio() #u*
	
	forcTr = forcTrio()
	
	Uerror = Uerr(Uth,Utr)
	
	#copie du dernier fichier Moyennes_spatiales_vitesse_rho_mu
	#recuperation du fichier
	#derniere ligne du ls
	ficLS = os.popen('ls -rt Moyennes_spatiales_vitesse_rho_mu_*')
	lignes = ficLS.readlines()
	Ligne = lignes[-1]
	#suppression du \n en fin de nom
	nomFic = Ligne[:len(Ligne)-1]
	os.popen('cp '+nomFic+' Moyennes_spatiales_vitesse_rho_mu')
	
	#ecriture des fichiers gnuplot et tableaux
	ecritureFichier(Uth,forcTh,Utr,forcTr,Uerror,nu)
