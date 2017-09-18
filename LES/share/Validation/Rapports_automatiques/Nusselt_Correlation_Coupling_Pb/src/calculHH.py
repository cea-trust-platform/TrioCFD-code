import os,os.path, sys, math
###
###
###
def properties(PG):
	# ouverture des fichiers
	fic = open(PG,'r')
	
	# lecture de ligne -> entetes
	fichier = fic.readlines()
	
	ligne = fichier [0]
	tLigne = ligne.split()
	mu = float(tLigne[0])
	cp = float(tLigne[3])
	lambada = float(tLigne[2])
	rho = float(tLigne[1])
	lx = float(tLigne[8])
	ly = float(tLigne[9])
	lz = float(tLigne[10])
	v = float(tLigne[11])
	nx = float(tLigne[12])
	nz = float(tLigne[13])
	lambda_c = float(tLigne[6])
	dhy = float(tLigne[14])
	fic.close()
	return mu,cp,lambada,rho,lx,ly,lz,v,nx,nz,lambda_c, dhy
##
#########
#########
def TS(ficTS):
	p = 0
#
	f1 = open(ficTS, 'r')
	lignes = f1.readlines()
	f1.close()
	ligne = lignes[1]
	ligne = ligne.strip()
	tLigne = ligne.split()
	valds = abs(float(tLigne[5]))
	ligne = lignes[-1]
	ligne = ligne.strip()
	tLigne = ligne.split()
	valTS = float(tLigne[1])

	return valds, valTS
#
def TF(ficTF):
	f2 = open(ficTF, 'r')
	lignes = f2.readlines()
	f2.close()
	ligne = lignes[-1]
	ligne = ligne.strip()
	tLigne = ligne.split()
	valTF = float(tLigne[1])

	return valTF
#
def PHI(ficPHI):
	f3 = open(ficPHI, 'r')
	fin = False
	ligne_dm = -1
	val_PHI = -1.
	while not fin:
		ligne = f3.readline()
		if ligne=='':
			fin = True
		else:
			ligne_dm += 1
			if ligne=='':
				fin = True
				ligne_dm = -1
			else:
				if (ligne.find('interface')>-1):
					ligne = f3.readline()
					ligneT = ligne.strip().split()
					val_PHI = abs(float(ligneT[-1]))
	f3.close()
	return val_PHI

def HH(valds, valTS, valTF, valPHI, lx, lz, nx, nz, lambda_c):
	Sface = lx/(nx-1) * lz/(nz-1)
	
	fi = open('h.dat','w')
	h = (valTS - valTF) * Sface / valPHI - valds / lambda_c
	h = 1 / h
	fi.write('%18.8f %18.8f\n' % (valTF, h))
	fi.close()
#
def Number(mu,cp,lambada,v,lx,lz,rho,dhy):
	Pr = mu*cp/lambada
	Re = v*dhy*rho/mu
	Nu = 0.023*(Re**0.8)*(Pr**(0.3333))
	hh = Nu*lambada/dhy
	fn = open('nh.dat','w')
	fn.write('%18.3f %18.f %18.2f %18.2f\n' % (Pr, Re, Nu, hh))
#########
######### 

if __name__=='__main__':
	#recuperation des parametres passes en ligne de commande
	args = sys.argv
	if len(args)!=5:
		print 'Erreur sur le nb d\'arguments fournis : Usage\npython f1 f2 f3 f4'
		sys.exit()
#
	ficTS = args[1]
	ficTF = args[2]
	ficPHI = args[3]
	PG = args[4]
#
	mu,cp,lambada,rho,lx,ly,lz,v,nx,nz,lambda_c,dhy = properties(PG)
	valds, valTS = TS(ficTS)
	valTF = TF(ficTF)
	valPHI = PHI(ficPHI)
	HH(valds, valTS, valTF, valPHI, lx, lz, nx, nz, lambda_c)
	Number(mu,cp,lambada,v,lx,lz,rho,dhy)
