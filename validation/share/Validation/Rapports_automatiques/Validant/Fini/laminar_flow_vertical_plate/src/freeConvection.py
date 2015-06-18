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
	TempLim=float(tLigne[10])
	rho=float(tLigne[8])
	
	fic.close()
	return mu,GrL,W,L,Temp,nbMailles,Cp,Pr,TempLim,rho
	
def ecritureFichier(nomFicTemp,nomFicVit,x,fluid,mu,GrL,L,Temp,TempLim,rho):
	#ecriture du fichier pour gnuplot
	nomFicGnuplot = 'freeConvection.dat'
	
	ficGnuPlot = open(nomFicGnuplot, 'w')
	ficTemp = open(nomFicTemp, 'r')
	ficVit = open(nomFicVit, 'r')

	ligneTemp = ficTemp.readline()
	ligneVit = ficVit.readline()
	
	deltaT = float(TempLim)-float(Temp)
	
	nu = float(mu)/float(rho)
	
	fin=False
	while not fin:
		ligneTemp = ficTemp.readline()
		ligneVit = ficVit.readline()
		if not ligneTemp or not ligneVit:
			fin=True
		else:
			tLigneTemp=ligneTemp.split()
			tLigneVit=ligneVit.split()
			y=tLigneTemp[0]
			T=tLigneTemp[1]
			u=tLigneVit[2]
			
			if fluid=="QC":
				T=(float(T)-float(Temp))/deltaT
					
			GrX = float(GrL)*math.pow((float(x)/float(L)),3)

			Y=(float(y)/float(x))*math.pow((GrX/4),(0.25))
			
			U=(float(u)*float(x))/(2*nu*math.pow(GrX,0.5))
			
			ficGnuPlot.write(' %18.4f %18.4f %18.4f\n' % (Y,float(T),U))
	
	ficGnuPlot.close()
	ficTemp.close()
	ficVit.close()

if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	mu,GrL,W,L,Temp,nbMailles,Cp,Pr,TempLim,rho = properties()
	
	#ecriture du fichier gnuplot
	ecritureFichier(args[0],args[1],args[2],args[3],mu,GrL,L,Temp,TempLim,rho)
