#Script de recuperation des proprietes physiques et de la geometrie

import os, sys, math

def readPropertiesData(nomFic):
	#initialisation
	properties = {}
	properties['mu'] = -1
	properties['Prandtl'] = -1
	properties['Cp'] = -1
	properties['gamma'] = -1
	properties['pression'] = -1	
	properties['D'] = -1
	properties['P'] = -1
	properties['L'] = -1
	properties['Tin'] = -1
	properties['Pvol'] = -1
	properties['Uin'] = -1
	properties['gg'] = -1
	properties['Q'] = -1
	# ouverture des fichiers
	fic = open(nomFic,'r')
	
#
	
	for ligne in fic:
		ligne = ligne.strip().lower()
		tLigne = ligne.split()
		if ligne.startswith('mu'):
			properties['mu'] = float(tLigne[-1])
		elif ligne.startswith('prandtl'):
			properties['Prandtl'] = float(tLigne[-1])
		elif ligne.startswith('cp'):
			properties['Cp'] = float(tLigne[-1])
		elif ligne.startswith('gamma'):
			properties['gamma'] = float(tLigne[-1])
		elif ligne.startswith('pression'):
			properties['pression'] = float(tLigne[-1])
		elif ligne.startswith('longueurs'):
			properties['D'] = float(tLigne[-3])
			properties['P'] = float(tLigne[-1])
			properties['L'] = float(tLigne[-2])
#
		if ligne.find('frontiere_ouverte_temperature_imposee')>-1:
			properties['Tin'] = float(tLigne[-1])
		elif ligne.find('puissance_thermique')>-1:
			properties['Pvol'] = float(tLigne[-2])
#
	properties['Uin'] = 3.72
	properties['gg'] = 9.81
	properties['Q'] = 0.00239
	fic.close()
	return properties


def ecritureFichier(properties):
	#ecriture du fichier
	nomFic = 'propertiesGeometry.dat'
	fichier = open(nomFic, 'w')
	fichier.write('%18.7f %18.3f %18.3f %18.3f %18.3f %18.5f %18.5f %18.5f %18.3f %18.3f %18.3f %18.3f %18.5f\n' % ( properties['mu'], properties['Prandtl'], properties['Cp'], properties['gamma'], properties['pression'], properties['D'], properties['P'], properties['L'], properties['Tin'], properties['Pvol'], properties['Uin'], properties['gg'], properties['Q']))
	fichier.close()

def getPropertiesFromdat():
	properties = {}
	nomFichier = 'propertiesGeometry.dat'
	if os.path.isfile(nomFichier):
		#recupere les donnees du fichier
		f = open(nomFichier, 'r')
		lignes = f.readlines()
		f.close()
	else:
		print 'Erreur getPropertiesFromdat : fichier %s non trouve !' % (nomFichier)
		sys.exit()
	ligne = (lignes[0]).strip()
	tabLigne = ligne.split()
	ind = 0
	try:
		properties['mu'] = float(tabLigne[ind])
		ind += 1
		properties['Prandtl'] = float(tabLigne[ind])
		ind += 1
		properties['Cp'] = float(tabLigne[ind])
		ind += 1
		properties['gamma'] = float(tabLigne[ind])
		ind += 1
		properties['pression'] = float(tabLigne[ind])
		ind += 1
		properties['D'] = float(tabLigne[ind])
		ind += 1
		properties['P'] = float(tabLigne[ind])
		ind += 1
		properties['L'] = float(tabLigne[ind])
		ind += 1
		properties['Tin'] = float(tabLigne[ind])
		ind += 1
		properties['Pvol'] = float(tabLigne[ind])
		ind += 1
		properties['Uin'] = float(tabLigne[ind])
		ind += 1
		properties['gg'] = float(tabLigne[ind])
		ind += 1
		properties['Q'] = float(tabLigne[ind])
	except IndexError:
		print 'Erreur getPropertiesFromdat : lecture element %d pour 0-%d elements...' % (ind, len(tabLigne)-1)
		sys.exit()
	except ValueError:
		print 'Erreur getPropertiesFromdat : lecture element %d n\'est pas un float (%s)...' % (ind, tabLigne[ind])
		sys.exit()
	return properties



if __name__ == '__main__':
	
	#recuperation du fichier data
	import glob
	#derniere ligne du ls
	#ficLS = os.popen('ls *.data')
	#lignes = ficLS.readlines()
	#Ligne = lignes[0]
	#suppression du \n en fin de nom
	#nomFic = Ligne[:len(Ligne)-1]
	listFics = glob.glob('*.data')
	if len(listFics)>0:
		nomFic = listFics[0]
		properties = readPropertiesData(nomFic)
	
		#ecriture du fichier
		ecritureFichier(properties)
	else:
		print 'Erreur propertiesGeometry : pas de fichier data trouve !'
