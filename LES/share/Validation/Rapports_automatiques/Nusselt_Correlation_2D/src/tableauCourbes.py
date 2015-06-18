#Script de construction des tableaux et des courbes

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
	H=float(tLigne[0])
	L=float(tLigne[1])
	mu=float(tLigne[2])
	lamb=float(tLigne[3])
	Cp=float(tLigne[4])	
	rho=float(tLigne[5])
	Tw=float(tLigne[6])
	D=float(tLigne[7])
	u=float(tLigne[8])
	
	fic.close()
	return H,L,mu,lamb,Cp,rho,Tw,D,u

#fonction de recuperation de T1 dans le cas VDF
def getT1():
	# ouverture des fichiers
	nomFic = 'LoiParoiNusseltImpose_VDF_COUPE_HE_T.coupe'
			
	fic = open(nomFic,'r')
	
	# lecture de la premiere ligne -> vide
	ligne = fic.readline()
	
	#recuperation de la deuxieme ligne
	ligne = fic.readline()
	
	tLigne = ligne.split()
	T1=float(tLigne[-1])
	fic.close()
	
	return T1

#fonction de recuperation de Flux dans le cas passe en parametre
def getFlux(discr):
	# ouverture des fichiers
	nomFic = 'LoiParoiNusseltImpose_'+discr+'_pbf_Diffusion_chaleur.face'
			
	fic = open(nomFic,'r')
	
	# lecture de la premiere ligne
	fichier = fic.readline()
	
	#on prend seulement le dernier temps (dernier bloc de donnees)
	commentaire=""
	i=-1
	while commentaire!="Flux":
		ligne = fic.readline()
		tLigne = ligne.split()
		commentaire = tLigne[1]
		i=i+1
	
	#lecture de l ensemble des donnees restantes
	tFic = fic.readlines()
	
	#on recupere la premiere ligne du dernier bloc (dernier temps)
	ligne = tFic[-i]
	tLigne=ligne.split()	
	Flux=float(tLigne[-1])
	
	fic.close()
	
	return Flux

#fonction de recuperation de T3 dans le cas VEF
def getT3():
	# ouverture des fichiers
	nomFic = 'LoiParoiNusseltImpose_VEF_COUPE_HE_T3.coupe'
			
	fic = open(nomFic,'r')
	
	# lecture de la premiere ligne -> vide
	ligne = fic.readline()
	
	#recuperation de la troisieme ligne
	ligne = fic.readline()
	ligne = fic.readline()
	
	tLigne = ligne.split()
	T3=float(tLigne[-1])
	fic.close()
	
	return T3

#fonction de recuperation de T2 dans le cas VEF
def getT2():
	# ouverture des fichiers
	nomFic = 'LoiParoiNusseltImpose_VEF_COUPE_HE_T2.coupe'
			
	fic = open(nomFic,'r')
	
	# lecture de la premiere ligne -> vide
	ligne = fic.readline()
	
	#recuperation de la deuxieme ligne
	ligne = fic.readline()
	
	tLigne = ligne.split()
	T2=float(tLigne[-1])
	fic.close()
	
	return T2

def ecritureFichier(hTheo,Flux,Tw,T1,hTrio,error):

	#ecriture du fichier pour la ligne du tableau
	nomFicLigne = 'ligneTableau.dat'
	fichier = open(nomFicLigne, 'w')
	fichier.write('%18.8f %18.8f %18.8f %18.8f %18.8f %18.2f\n' % (hTheo,Flux,Tw,T1,hTrio,error))
	fichier.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	discr=args[0]
	
	#recuperation des donnees
	H,L,mu,lamb,Cp,rho,Tw,D,u = properties()
	
	#theo
	nu = float(mu)/float(rho)

	Pr = nu*float(rho)*float(Cp)/float(lamb)

	Re = float(u)*float(D)/nu

	Nusselt = 0.023*math.pow(Re,0.8)*math.pow(Pr,1.0/3.0)

	hTheo = Nusselt*float(lamb)/float(D)

	#trio	
	S=0.001
	
	if discr=="VDF":
		T1 = getT1()
	else :
		T2 = getT2()
		T3 = getT3()
		T1 = 0.5*(float(T2)+float(T3))

	Flux = getFlux(discr)
	
	hTrio=float(Flux)/((float(Tw)-float(T1))*float(S))
	
	error = abs((hTheo - hTrio) / hTheo)*100
	
	#ecriture des fichiers gnuplot et tableaux
	ecritureFichier(hTheo,Flux,Tw,T1,hTrio,error)
