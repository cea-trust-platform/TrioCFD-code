#Script de calcul de Y,u et T*

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
	GrL=tLigne[1]
	W=float(tLigne[2])
	L=float(tLigne[3])
	Temp=float(tLigne[4])
	nbMailles=float(tLigne[5])
	Cp=float(tLigne[6])
	Pr=float(tLigne[7])
	rho=float(tLigne[8])
	U=float(tLigne[9])
	TempLim=float(tLigne[10])
	
	
	fic.close()
	return mu,GrL,W,L,Temp,nbMailles,Cp, Pr,rho,U,TempLim

def ecritureFichier(L,nbMailles,Cp,fluid,GrL,Pr,case,mu,rho,U,Temp,TempLim):
	#ecriture du fichier pour gnuplot
	nomFicGnuplot = 'mixedConvection.dat'
	nomFicRead = 'MixConv_pb_Diffusion_chaleur.face'
	if fluid=="QC" and case=="VDF": 
		nomFicRead = 'Mixconv_QC_pb_Diffusion_chaleur.face'
	
	ficGnuPlot = open(nomFicGnuplot, 'w')
	ficRead = open(nomFicRead, 'r')

	#lecture de la premiere ligne
	ligne = ficRead.readline()
	
	#on prend seulement le dernier temps (dernier bloc de donnees)
	commentaire=""
	i=-1
	while commentaire!="Flux":
		ligne = ficRead.readline()
		tLigne = ligne.split()
		commentaire = tLigne[1]
		i=i+1
	
	tFic = ficRead.readlines()
	
	#calcul de lambda
	lamb = float(mu)*float(Cp)/float(Pr)
	
	#calcul de la vitesse ReL
	ReL = float(rho)*float(U)*float(L)/float(mu)

	#calcul du deltaT
	deltaT=float(TempLim)-float(Temp)
	
	while i>0:
			ligne = tFic[-i]
			tLigne=ligne.split()
			x=float(tLigne[6])
			w=float(tLigne[-1])
			
			X=float(x)/float(L)
			
			Nux = w*float(L)/(deltaT/float(nbMailles))*(float(x)/float(lamb))
			
			GrX = (float(GrL)*math.pow(float(x),3))/math.pow(float(L),3)
			
			Rex = float(ReL)*float(x)/float(L)
			
			Rix = GrX/(Rex*Rex)
			
			NuxRex = Nux/math.pow(Rex,0.5)
			
			ficGnuPlot.write('%18.4f %18.4f %18.4f %18.4f\n' % (X,Nux,Rix,NuxRex))
			
			i=i-1
			
	ficGnuPlot.close()
	ficRead.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	mu,GrL,W,L,Temp,nbMailles,Cp, Pr,rho,U,TempLim = properties()
	
	#ecriture du fichier gnuplot
	ecritureFichier(L,nbMailles,Cp,args[0],GrL,Pr,args[1],mu,rho,U,Temp,TempLim)
