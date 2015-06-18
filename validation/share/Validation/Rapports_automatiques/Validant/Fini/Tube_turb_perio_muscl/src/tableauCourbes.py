#Script de construction des tableaux et des courbes
import os

import optparse
import math

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
	U=float(tLigne[2])
	Z=10	
	L=9
	D=4
	fic.close()
	return D,mu,rho,U,Z,L

#fonction de recuperation de la colonne y+
def yPlus():
	# ouverture des fichiers
	nomFicUstar = 'TUBE_TURB_perio_pb_Ustar.face'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	#on enleve les 5 dernieres lignes avant le tableau
	i=0
	while i < 5: 
		del fichier [-1]
		i=i+1
	
	commentaire=""
	i=len(fichier)
	cumul=0
	y=0.
	
	valz=float(L/2.)
	#on va prendre la valeur des moyennes de y+ entre L/2 et L (eviter l'entree)
	while commentaire!="----------------":
		ligne=fichier[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la troisieme colonne (z)
				z=float(tLigne[2])
				if z > valz :
					y = y + float(tLigne[4])
					cumul = cumul + 1
		i=i-1			
	
	y=y/float(cumul)
	ficUstar.close()
	
	return y

def uPlus():
	# ouverture des fichiers
	nomFicUstar = 'TUBE_TURB_perio_pb_Ustar.face'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	#on enleve les 5 dernieres lignes avant le tableau
	i=0
	while i < 5: 
		del fichier [-1]
		i=i+1
	
	commentaire=""
	i=len(fichier)
	cumul=0
	u=0.
	
	valz=float(L/2.)
	#on va prendre la valeur des moyennes de y+ entre L/2 et L (eviter l'entree)
	while commentaire!="----------------":
		ligne=fichier[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la troisieme colonne (z)
				z=float(tLigne[2])
				if z > valz :
					u = u + float(tLigne[3])
					cumul = cumul + 1
		i=i-1			
	
	u=u/float(cumul)
	ficUstar.close()
	
	return u

def Utrio():
	# ouverture des fichiers
	nomFicUstar = 'TUBE_TURB_perio_pb_Ustar.face'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	#on enleve les 5 dernieres lignes avant le tableau
	i=0
	while i < 5: 
		del fichier [-1]
		i=i+1
	
	commentaire=""
	i=len(fichier)
	cumul=0
	utau=0.
	
	valz=float(L/2.)
	#on va prendre la valeur des moyennes de y+ entre L/2 et L (eviter l'entree)
	while commentaire!="----------------":
		ligne=fichier[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la troisieme colonne (z)
				z=float(tLigne[2])
				if z > valz :
					utau = utau + float(tLigne[5])
					cumul = cumul + 1
		i=i-1			
	
	utau=utau/float(cumul)
	ficUstar.close()
	
	return utau
	
#fonction permettant de calculer la vitesse moyenne calculee par le code (differente de Ucible?)
def UbulkCalc(D):
	#calcul de la la surface 
	Surf=3.14159*float(D)*float(D)/4.
	#ouverture des fichiers
	nomFicUbulk = 'TUBE_TURB_perio_pb_Debit.out'
        ficUbulk = open(nomFicUbulk,'r')

        fichier = ficUbulk.readlines()
        #derniere ligne
        ligne = fichier[-1]
        tLigne = ligne.split()
        Ucalc=tLigne[3]
        Ucalc=float(Ucalc)/(float(Surf))
        ficUbulk.close()
        return Ucalc

def Re(mu,rho,D,U):
	nu = float(mu)/float(rho)
	reh = float(U)*float(D)/nu
	return reh

def Utheo(U,reh):
	lamb = 0.316*math.pow(reh,-0.25)
	utau = float(U)*math.sqrt(float(lamb)/8.00)
	return utau

def Uerr(Utheo,Utrio):
	Uerror = float(abs(Utheo-Utrio)/Utheo)*100
	return Uerror

def DHtheo(Utheo,D):
	DHth = (4*math.pow(float(Utheo),2))/float(D)
	return DHth

def DHtrio(Utrio,D):
	DHtrio=math.pow(Utrio,2)*4/float(D)
	return DHtrio


def DHtrio2():
	# ouverture des fichiers
	nomFic = os.popen('ls -rt *Pressure_Gradient_pb_perio').readlines()[-1]
	nomFic = nomFic[:len(nomFic)-1] # Suppress /n
			
	fic = open(nomFic,'r')
	
	# lecture de ligne -> entetes
	fichier = fic.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	DHtrio=float(tLigne[1])
	fic.close()
		
	return DHtrio
		
def DHerr(DeltaHTheo,DeltaHTrio):
	DeltaHError = float(abs(DeltaHTheo-DeltaHTrio)/DeltaHTheo)*100
	return DeltaHError

def ecritureFichier(D,mu,rho,U,Ucalc,reh,rehCalc,yPlus,uPlus,Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2,DRey):

	nu = float(mu)/float(rho)

	#ecriture du fichier pour les vitesses et reynolds
	nomFicRe = 'ligneReynolds.dat'
	fichier = open(nomFicRe, 'w')
	fichier.write('%18.8f %18.8f %18.1f %18.1f %18.2f\n' % (U,Ucalc,reh,rehCalc,DRey))
	fichier.close()

#ecriture du fichier pour la ligne du tableau
	nomFicLigne = 'ligneTableau.dat'
	fichier = open(nomFicLigne, 'w')
	fichier.write('%18.8f %18.8f %18.2f %18.8f %18.8f %18.2f %18.8f %18.2f\n' % (Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2))
	fichier.close()

	#ecriture du fichier pour le premier point
	nomFicFirstPoint = 'first_Point.dat'
	fichier = open(nomFicFirstPoint, 'w')
	fichier.write('%18.8f %18.8f\n' % (yPlus,uPlus))
	fichier.close()
	 
	#ecriture du fichier pour la courbe Axial velocity	
	nomFicRead = 'TUBE_TURB_perio_SONDE_VP.coupe'
	nomFic = 'axial_velocity.dat'
	ficRead = open(nomFicRead, 'r')
	#1ere ligne vide
	ligne=ficRead.readline()
	fichier = open(nomFic, 'w')
	fin = False
	while not fin:
		ligne = ficRead.readline()
		if not ligne:
			fin=True
		else:
			tLigne = ligne.split()
			y = tLigne[0]
			u = tLigne[-1]
			yAdim = (float(y)+2)*float(Utr)/nu
			uAdim = float(u)/float(Utr)
			fichier.write('%18.8f %18.8f\n' % (yAdim,uAdim))
	fichier.close()
	ficRead.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	#recuperation des donnees
	D,mu,rho,U,Z,L = properties()

	reh = Re(mu,rho,D,U)
#	Ucalc = 27.5
	Ucalc = UbulkCalc(D)
	rehCalc = Re(mu,rho,D,Ucalc)
	y = yPlus()
	u = uPlus()
	Uth = Utheo(U,reh)
	Uthcalc = Utheo(Ucalc,rehCalc)
	Utr = Utrio() #u*
	Uerror = Uerr(Uth,Utr)
	DeltaHTheo = DHtheo(Uth,D)
	DeltaHTrio = DHtrio(Utr,D)
	DeltaHTrio2 = DHtrio2()
	DeltaHError = DHerr(DeltaHTheo,DeltaHTrio)
	DeltaHError2 = DHerr(DeltaHTheo,DeltaHTrio2)
	DRey = DHerr(reh,rehCalc)
	
	#ecriture des fichiers gnuplot et tableaux
	ecritureFichier(D,mu,rho,U,Ucalc,reh,rehCalc,y,u,Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2,DRey)
	
