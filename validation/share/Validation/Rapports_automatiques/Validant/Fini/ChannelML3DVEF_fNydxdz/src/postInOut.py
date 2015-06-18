# Code pour le post-traitement du cas unique en Entree Sortie de la fiche

import optparse
import os
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
	h=float(tLigne[2])
	L=float(tLigne[3])
	Ub=float(tLigne[4])
	
	fic.close()
	return mu,rho,h,L,Ub


def yPlus(L):
	# ouverture des fichiers
	nomFicUstar = '3D_keps_pb_Ustar.face'
			
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
	
	valx=float(0.75*L)
	#on va prendre la valeur des moyennes de y+ entre 3/4 L et L (eviter l'entree)
	while commentaire!="----------------":
		ligne=fichier[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la premiere colonne (x)
				x=float(tLigne[0])
				if x > valx :
					y = y + float(tLigne[4])
					cumul = cumul + 1
		i=i-1			
	
	y=y/float(cumul)
	ficUstar.close()
	
	return y

def u_tau(L):
	# ouverture des fichiers
	nomFicUstar = '3D_keps_pb_Ustar.face'
			
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
	u_tau=0.
	
	valx=float(0.75*L)
	#on va prendre la valeur des moyennes de y+ entre 3/4L et L (eviter l'entree)
	while commentaire!="----------------":
		ligne=fichier[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la premiere colonne (x)
				x=float(tLigne[0])
				if x > valx :
					u_tau = u_tau+ float(tLigne[5])
					cumul = cumul + 1
		i=i-1			
	
	u_tau=u_tau/float(cumul)
	ficUstar.close()
	
	return u_tau

def retau(rho,h,mu,u_tau):

	Re=rho*10*h/mu/2
	Puissance = math.pow(float(Re),(7.0/8.0))
	Re_tau_Theo=0.1754*Puissance

	Re_tau=rho*u_tau*(h/2)/mu	

	return Re_tau,Re_tau_Theo


def perte_charge(mu,rho,h,Ub,L):
	# ouverture des fichiers
	nomFic = os.popen('ls -rt *Pressure_Gradient_pb_periox').readlines()[-1]
	nomFic = nomFic[:len(nomFic)-1] # Suppress /n
		
	fic = open(nomFic,'r')
	
	# lecture de ligne -> entetes
	fichier = fic.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	Perte=float(tLigne[1])*float(L)
	fic.close()
	Re=rho*Ub*(h/2)/mu
	Puissance = math.pow(float(Re),(7.0/8.0))
	Pertetheo=(0.1754*Puissance*mu)*(0.1754*Puissance*mu)
	
	return Perte,Pertetheo



def ecritureFichier(y_plus,u_tau_trio,Re_tau,Re_tau_theo,ReError):

	nomFic = 'ligneTableau.dat'
	fichier = open(nomFic, 'w')
	fichier.write('%18.4f %18.4f %18.4f %18.4f %18.4f\n' % (y_plus,u_tau_trio,Re_tau,Re_tau_theo,ReError)) 
	fichier.close()

def ecrire_profils_vitesse(sonde,u_tau,mu,rho):

	ficRead = open(sonde, 'r')
	
	nomFic = 'Velocity_profile.dat'
	ficWrite = open(nomFic, 'w')

	nomFic2 = 'Dimensionless_velocity_profile.dat'
	ficWrite2 = open(nomFic2, 'w')
	
	ligne = ficRead.readlines()	
	n=len(ligne)
	
	for i in range (1,n):
			
		tligne=ligne[i].split()			
		ficWrite.write('%18.8f %18.8f\n' % (float(tligne[0]),float(tligne[1])))
		ficWrite2.write('%18.8f %18.8f\n' % (float(tligne[0])*u_tau_trio*rho/mu,float(tligne[1])/u_tau_trio))
	
	ficWrite.close()
	ficWrite2.close()
	ficRead.close()

def ecrire_profils_visco(sonde,u_tau,h):

	ficRead = open(sonde, 'r')
	
	nomFic = 'Dimensionless_turbulent_viscosity.dat'
	ficWrite = open(nomFic, 'w')
	
	ligne = ficRead.readlines()	
	n=len(ligne)
	
	for i in range (1,n):
			
		tligne=ligne[i].split()	
		
		if float(tligne[0])<1.0:		
			ficWrite.write('%18.8f %18.8f\n' % (float(tligne[0])/(0.5*h),float(tligne[1])/(u_tau*(h/2))))
		
	
	ficWrite.close()

	ficRead.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	mu,rho,h,L,Ub = properties()
	
	#Calcul de Retau et ReTautheo
	
	y_plus = yPlus(L)
	u_tau_trio=u_tau(L)
	Re_tau,Re_tau_theo=retau(rho,h,mu,u_tau_trio)
	ReError=(abs(Re_tau-Re_tau_theo)/Re_tau_theo)*100
		
	#ecriture du fichier ligneTableau.dat
	ecritureFichier(y_plus,u_tau_trio,Re_tau,Re_tau_theo,ReError)

	#ecriture du profils de vitesse Vx=f(y) et U+=f(y+) en x=39
	fichierSonde='3D_keps_SONDE_VIT.coupe'
	ecrire_profils_vitesse(fichierSonde,u_tau_trio,mu,rho)


	#ecriture du profils de viscosite turbulente adimensionnee
	fichierSonde2='3D_keps_SONDE_VISC_TURB.coupe'
	ecrire_profils_visco(fichierSonde2,u_tau_trio,h)

















