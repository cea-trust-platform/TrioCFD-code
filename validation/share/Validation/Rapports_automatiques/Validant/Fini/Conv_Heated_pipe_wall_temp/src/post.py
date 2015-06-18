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
	D=float(tLigne[0])
	mu=float(tLigne[1])
	lamb=float(tLigne[2])
	Cp=float(tLigne[3])
	rho=float(tLigne[4])	
	U=float(tLigne[5])
	Qw=float(tLigne[6])
	Z=float(tLigne[7])
	
	fic.close()
	return D,mu,lamb,Cp,rho,U,Qw,Z

#fonction de recuperation de la colonne y+
def yPlus():
	# ouverture des fichiers
	nomFicUstar = 'Conduite_pb_Ustar.face'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#en periodique on prend la valeur moyenne calculee en fin de fichier
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	y=float(tLigne[3])
	ficUstar.close()
	
	return y

def uPlus():
	# ouverture des fichiers
	nomFicUstar = 'Conduite_pb_Ustar.face'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	u=float(tLigne[1])
	ficUstar.close()
	
	return u

def Utrio():
	# ouverture des fichiers
	nomFicUstar = 'Conduite_pb_Ustar.face'
			
	ficUstar = open(nomFicUstar,'r')
	
	# lecture de ligne -> entetes
	fichier = ficUstar.readlines()
	
	#tant que la derniere ligne est vide
	while fichier[-1]=="" or fichier[-1]=="\n":
		del fichier [-1]
	
	ligne = fichier [-1]
	tLigne = ligne.split()
	utau=float(tLigne[5])
	ficUstar.close()
	
	return utau

#fonction permettant de calculer la vitesse moyenne calculee par le code (differente de Ucible?)
def UbulkCalc(rho,D):
	#calcul de la la surface : attention dans RhoU_sortie surface double 
	Surf=3.14159*float(D)*float(D)/4.
	#ouverture des fichiers
	nomFicUbulk = 'RhoU_sortie'
        ficUbulk = open(nomFicUbulk,'r')

        fichier = ficUbulk.readlines()
        #derniere ligne
        ligne = fichier[-1]
        tLigne = ligne.split()
        Ucalc=tLigne[1]
        Ucalc=float(Ucalc)/float(rho)/(2.*float(Surf))
        ficUbulk.close()
        return Ucalc

	
def Tparoi(Z,col):
	# recherche de Tparoi au meme endroit que la sonde pour tracer les profils
	# ouverture des fichiers
	nomFicTw = 'Conduite_pb_Nusselt.face'
			
	ficTw = open(nomFicTw,'r')
	
	tFic = ficTw.readlines()
	commentaire=""
	i=len(tFic)
#	print "*********************** ",i
	
	valz=Z
	deltaz=Z
	#determination de la valeur la plus proche de Z (par valeur inf)
	while commentaire!="----------------":
		ligne=tFic[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la troisieme colonne (z)
				z=float(tLigne[2])
				if (Z-z) < abs(deltaz) and (Z-z)>0 :
					valz = z
					deltaz = Z-z
		i=i-1			
	#on cherche la valeur x la plus proche de la paroi (X=-1) pour Z 
	#attention on predefinit X=-1 (lire dans sonde, ce serait mieux)
	commentaire=""
	i=len(tFic)
	X=-1.
	valx=X
	deltaX=X
	while commentaire!="----------------":
		ligne=tFic[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				if float(tLigne[2])==valz :
					x = float(tLigne[0])
					if abs(X-x) <= abs(deltaX) :
						valx = x
						deltaX = X-x
		i=i-1
			
	#moyenne des valeurs Tparoi sur z = Z et x = X			
	commentaire=""
	i=len(tFic)
	cumul=0
	nbTemps=0
	while commentaire!="----------------":
		ligne=tFic[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				if float(tLigne[2])==valz and float(tLigne[0])==valx:
					# on lit la derniere valeur (Tp equiv)
					T=float(tLigne[col])
					cumul = cumul + T
					nbTemps = nbTemps + 1
		i=i-1
			
	Tw = cumul/float(nbTemps)
	
	ficTw.close()
	
	return Tw




def Tparoi_max_calcule(Z,col):
	# recherche de Tparoi au meme endroit que la sonde pour tracer les profils
	# ouverture des fichiers
	nomFicTw = 'Conduite_pb_Nusselt.face'
			
	ficTw = open(nomFicTw,'r')

	tFic = ficTw.readlines()
	commentaire=""
	i=len(tFic)

	valmax=0

	
	#moyenne des valeurs Tparoi sur z = Z et x = X			
	commentaire=""
	i=len(tFic)
	while commentaire!="----------------":
		ligne=tFic[i-1]
		if ligne != "\n" and ligne != "":
		        #on eclate la ligne en un tableau de valeurs
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")

				val=float(tLigne[col])
				if val > valmax :
					valmax=val

		i=i-1
			
	
	ficTw.close()
	
	return valmax



def temperature_physique(nomFicRead):
	
	ficRead = open(nomFicRead,'r')

	# On passe les 4 premieres lignes de commentaires
	ligne = ficRead.readline()
	ligne = ficRead.readline()
	ligne = ficRead.readline()
	ligne = ficRead.readline()

	# balayage du ficher pour ne garder que la valeur de T de la derniere ligne 
	fin = False
	while not fin:
		ligne = ficRead.readline()
		if not ligne:
			fin=True
		else:
			tLigne = ligne.split()
			t = tLigne[-1]
	ficRead.close()

	# On affecte la derniere valeur trouvee
	T_phys=float(t)
	return T_phys
	










	
def Ttotal():
	#en periodique on considere qu'on peut faire la moyenne sur toute la longueur
	# ouverture des fichiers
	nomFicTw = 'Conduite_pb_Nusselt.face'
			
	ficTw = open(nomFicTw,'r')

	
	
	tFic = ficTw.readlines()
	commentaire=""
	i=len(tFic)
	cumul=0
	nbTemps=0
	while commentaire!="----------------":
		ligne=tFic[i-1]
		if ligne != "\n" and ligne != "":
			tLigne = ligne.split("|")
			commentaire = tLigne[0]
			if commentaire!="----------------":
				#on eclate la ligne en un tableau de valeurs
				tLigne=ligne.split("|")	
				# on lit la derniere valeur
				T=float(tLigne[-1])
				cumul = cumul + T
				nbTemps = nbTemps + 1
		i=i-1
	
	Tw = cumul/float(nbTemps)
	
	ficTw.close()
	
	return Tw
	
def Tbulk():
	#ouverture des fichiers
	nomFicTbulk = 'Tmoyen_sortie'
			
	ficTbulk = open(nomFicTbulk,'r')
	
	fichier = ficTbulk.readlines()
	#derniere ligne
	ligne = fichier[-1]
	tLigne = ligne.split()
	
	Tbulk=tLigne[1]
	
	ficTbulk.close()
	
	return Tbulk


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
	nomFic = os.popen('ls -rt *Pressure_Gradient_pb_sortie').readlines()[-1]
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

def ecritureFichier(D,mu,lamb,Cp,rho,U,Ucalc,reh,rehCalc,yPlus,uPlus,Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2):

	nu = float(mu)/float(rho)

	#ecriture du fichier pour les vitesses et reynolds
	nomFicRe = 'ligneReynolds.dat'
	fichier = open(nomFicRe, 'w')
	fichier.write('%18.8f %18.8f %18.8f %18.2f \n' % (U,Ucalc,reh,rehCalc))
	fichier.close()

	#ecriture du fichier pour la ligne du tableau
	nomFicLigne = 'ligneTableau.dat'
	fichier = open(nomFicLigne, 'w')
	fichier.write('%18.8f %18.8f %18.2f %18.8f %18.8f %18.2f %18.8f %18.2f %18.2f\n' % (Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2,yPlus))
	fichier.close()

	#ecriture du fichier pour le premier point
	nomFicFirstPoint = 'first_Point.dat'
	fichier = open(nomFicFirstPoint, 'w')
	fichier.write('%18.8f %18.8f\n' % (yPlus,uPlus))
	fichier.close()
	 
	#ecriture du fichier pour la courbe Axial velocity	
	nomFicRead = 'Conduite_SONDE_VP.coupe'
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
			yAdim = (float(y)+1)*float(Utr)/nu
			uAdim = float(u)/float(Utr)
			fichier.write('%18.8f %18.8f\n' % (yAdim,uAdim))
	fichier.close()
	ficRead.close()

def temperatureProfile(Tw,Qw,mu,rho,Cp,Utr):
	
	nu = float(mu)/float(rho)
	nomFicRead = 'Conduite_SONDE_TP.coupe'
	nomFic = 'temperature_profile.dat'
	ficRead = open(nomFicRead,'r')
	fichier = open(nomFic,'w')
	#1ere ligne vide
	ligne=ficRead.readline()
	fin = False
	while not fin:
		ligne = ficRead.readline()
		if not ligne:
			fin=True
		else:
			tLigne = ligne.split()
			y = tLigne[0]
			t = tLigne[-1]
			yAdim = (float(y)+1)*float(Utr)/nu
			#yAdim = (float(y))*float(Utr)/nu
			#yAdim = (2-float(y))*float(Utr)/nu
			tAdim = (Tw-float(t))/(float(Qw)/(float(rho)*float(Cp)*float(Utr)))
			fichier.write('%18.8f %18.8f\n' % (yAdim,tAdim))
	fichier.close()
	ficRead.close()

	
def temperatureFirstPoint(Tfp,Tw,Qw,mu,rho,Cp,Utr,yPlus):
	
        #ecriture du fichier pour le premier point
	nomFic = 'temperature_first_point.dat'
	fichier = open(nomFic,'w')

	Ttau=Qw/(rho*Cp*Utr)
	Tplus=(Tw-Tfp)/Ttau
	
	fichier.write('%18.8f %18.8f\n' % (yPlus,Tplus))
	fichier.close()
	

	
def temperatures_paroi(Tfp,Tw,Tface,Tmax_paroi_post):
	
        #ecriture du fichier
	nomFic = 'temperatures_paroi.dat'
	fichier = open(nomFic,'w')

	#Tface: T sur la face de paroi
	#Tw  : T paroi equivalente
	#Tmax_paroi_post : T paroi post traitee max
        ecart_Tcalc=100.*abs(Tw-Tface)/Tw
	ecart_Tpost=100.*abs(Tw-Tmax_paroi_post)/Tw

	
	fichier.write('%14.6f %14.2f %14.2f %14.2f %14.2f %14.2f\n' % (Tfp,Tw,Tface,Tmax_paroi_post,ecart_Tcalc,ecart_Tpost))
	fichier.close()



def kader(mu,Cp,lamb,rho,Pr,Utr):
	
	nu = float(mu)/float(rho)
	
	nomFic = 'kader.dat'
	fichier = open(nomFic,'w')
	
	y=5
	while y<=10000:
		
		beta = math.pow(3.85*math.pow(Pr,1.0/3.0)-1.3,2)+2.12*math.log(Pr)
		gamma = 0.01*math.pow(Pr*y,4)/(1+5*y*math.pow(Pr,3))
		tAdim = Pr*y*math.exp(-1*gamma)+(2.12*math.log(1+y)+beta)*math.exp(-1.0/gamma)
		fichier.write('%18.8f %18.8f\n' % (y,tAdim))
		y=y+10

		
	fichier.close()





def PryPlus(Pr,Utr):
	nu = float(mu)/float(rho)
	
	nomFic = 'pryPlus.dat'
	fichier = open(nomFic,'w')
	
	y=1
	while y<=100:
		
		y=y+1
		tAdim = Pr*y
		fichier.write('%18.8f %18.8f\n' % (y,tAdim))
	
	fichier.close()

def nusselt(TwTotal,Qw,D,lamb,Tbulk,reh,Pr,retr,Tface,Tmax_paroi_post):
	nomFic = 'nusselt.dat'
	fichier = open(nomFic,'w')
	nomFic2 = 'calcul_nusselt.dat'
	fichier2 = open(nomFic2,'w')
	nomFic3 = 'nusselt_T_post.dat'
	fichier3 = open(nomFic3,'w')

	
	Nu_trio=float(float(Qw)*float(D)/float(lamb)/(float(TwTotal)-float(Tbulk)))
	Nu_trio_face=float(float(Qw)*float(D)/float(lamb)/(float(Tface)-float(Tbulk)))
	Nu_trio_post=float(float(Qw)*float(D)/float(lamb)/(float(Tmax_paroi_post)-float(Tbulk)))

	Nu_Colburn = 0.023*math.pow(float(reh),0.8 )*math.pow(float(Pr),0.3333)
	Nu_Dittus  = 0.026*math.pow(float(reh),0.8 )*math.pow(float(Pr),0.3)

        var_col=float(abs(Nu_Colburn-float(Nu_trio))/Nu_Colburn)*100
        var_col_face=float(abs(Nu_Colburn-float(Nu_trio_face))/Nu_Colburn)*100
        var_col_post=float(abs(Nu_Colburn-float(Nu_trio_post))/Nu_Colburn)*100
        var_dit=float(abs(Nu_Dittus-float(Nu_trio))/Nu_Dittus)*100
	diffT=float(TwTotal)-float(Tbulk)

		
	fichier.write('%14.5f %14.2f %14.2f %12.2f %12.2f %12.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n' % (Pr, reh, retr, Nu_trio,Nu_Colburn,Nu_Dittus,var_col,var_dit,TwTotal,Qw,D,lamb,float(Tbulk)))
	fichier.close()
	fichier2.write('%10.2f %10.2f %10.4f %10.2f %10.2f %12.2f %12.2f\n' % (float(Qw),float(D),float(lamb),float(TwTotal),float(Tbulk),diffT,float(Nu_trio)))
	fichier2.close()
	fichier3.write('%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f \n' % (Nu_trio,Nu_trio_face,Nu_trio_post,Nu_Colburn,var_col,var_col_face,var_col_post))
	fichier3.close()




if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	(options, args) = parser.parse_args()
	
	# recuperation des donnees
	D,mu,lamb,Cp,rho,U,Qw,Z = properties()



	# Calcul de variables intermediaires
	Pr = float(mu)*float(Cp)/float(lamb)	
	reh = Re(mu,rho,D,U)
	Ucalc = UbulkCalc(rho,D)
	rehCalc = Re(mu,rho,D,Ucalc)
	
	y = yPlus()
	u = uPlus()
	
	Uth = Utheo(U,reh)
	Uthcalc = Utheo(Ucalc,rehCalc)
	Utr = Utrio() #u*
        retr = Re(mu,rho,D,Utr)
	Uerror = Uerr(Uth,Utr)
	DeltaHTheo = DHtheo(Uth,D)
	DeltaHTrio = DHtrio(Utr,D)
	DeltaHTrio2 = DHtrio2()
	DeltaHError = DHerr(DeltaHTheo,DeltaHTrio)
	DeltaHError2 = DHerr(DeltaHTheo,DeltaHTrio2)
	
	# temperature de paroi
	col=-1
	Tw = Tparoi(Z,col)
	Tw_max = Tparoi_max_calcule(Z,col)
	# temperature de la premiere maille de la sonde
	col=6
	Tfp = Tparoi(Z,col)
	Tfp_max = Tparoi_max_calcule(Z,col)
	# temperature de la face de la maille
	col=7
	Tface = Tparoi(Z,col)
	Tface_max = Tparoi_max_calcule(Z,col)

	# temperature moyenne equivalente sur la paroi
	TwTotal = Ttotal()
	# temperature bulk
	Tbulk = Tbulk()
	# temperature_physique maximum et moyenne 

	Fic_Tmax = 'Conduite_SONDE_T_PHYS_PAROI_MAX.son'
	T_phys_max=temperature_physique(Fic_Tmax)
	#print T_phys_max

	
	Fic_Tmoy = 'Conduite_SONDE_T_PHYS_PAROI_MOY.son'
	T_phys_moy=temperature_physique(Fic_Tmoy)
	#print T_phys_moy



	
	
	#ecriture des fichiers gnuplot et tableaux pour l'hydraulique
	ecritureFichier(D,mu,lamb,Cp,rho,U,Ucalc,reh,rehCalc,y,u,Uth,Utr,Uerror,DeltaHTheo,DeltaHTrio,DeltaHError,DeltaHTrio2,DeltaHError2)
	
	#ecriture du fichier gnuplot pour T+=f(y+)
	temperatureProfile(Tw,Qw,mu,rho,Cp,Utr)
	
	#ecriture du fichier gnuplot pour le trace de la loi de Kader
	kader(mu,Cp,lamb,rho,Pr,Utr)

	#ecriture du fichier donnant le premier point pour la temperature
	temperatureFirstPoint(Tfp,Tw,Qw,mu,rho,Cp,Utr,y)	
		
	#ecriture du fichier gnuplot pour la courbe Pry+
	PryPlus(Pr,Utr)
	
	#ecriture du nusselt
	#  Calcul du Nusselt avec Tequi
	#  Comparaison du Nusselt calcule avec Tequi avec Dittus et Colburn
	#  Comparaison des Nusselt calcules avec Tequi Tpost et Tface et Colburn

	nusselt(TwTotal,Qw,D,lamb,Tbulk,reh,Pr,retr,Tface,T_phys_moy)
	#nusselt(Tw_max,Qw,D,lamb,Tbulk,reh,Pr,retr,Tface_max,T_phys_max)
	
	
	#temperatures_paroi
	# Pour tableau de comparaison Tequi Tface et Tpost
	
	temperatures_paroi(Tfp,TwTotal,Tface,T_phys_moy)
	#temperatures_paroi(Tfp,Tw_max,Tface_max,T_phys_max)
	
