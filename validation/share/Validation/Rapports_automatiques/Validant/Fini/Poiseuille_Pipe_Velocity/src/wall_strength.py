import os,os.path, sys, math


def properties():
	# ouverture des fichiers
	nomFp = 'propertiesGeometry.dat'
			
	fic = open(nomFp,'r')
	
	# lecture de ligne -> entetes
	fichier = fic.readlines()
	
	ligne = fichier [0]
	tLigne = ligne.split()
	rho     = float(tLigne[0])
	mu      = float(tLigne[1])
	rayon   = float(tLigne[2])
	longueur= float(tLigne[3])
	vitesse = float(tLigne[4])
	fic.close()
	return rho,mu,rayon,longueur,vitesse
##
##


def getNormeF(nomFic, rho,mu,rayon,longueur,vitesse):
#
	pi = 3.141592654
	DP = -8*mu*vitesse / (rayon*rayon)
	f = open(nomFic, 'r')
	lignes = f.readlines()
	f.close()
	#for ligne in lignes:
	ligne = lignes[-1]
	ligne = ligne.strip()
	if ligne[0]!='#':
		tLigne = ligne.split()
		valtemps = float(tLigne[0])
		fic = open('calculF.dat', 'w')
		fz = float(tLigne[3])
		uz = math.sqrt(fz / (rho*longueur*2*pi))
		uth = math.sqrt(-DP/(rho*2)) 
		Du = 100*abs(uz-uth)/max(uz,uth)
		fic.write('%18.3f %18.3f %18.2f\n' % (uth,uz,Du))
#
#
if __name__=='__main__':
	#recuperation des parametres passes en ligne de commande
	args = sys.argv
	if len(args)!=2:
		print 'Erreur sur le nb d\'arguments fournis : Usage  python wall_strength nomFichier.out'
		sys.exit()

	nomFic = args[1]
	rho,mu,rayon,longueur,vitesse = properties()
	getNormeF(nomFic, rho,mu,rayon,longueur,vitesse)



