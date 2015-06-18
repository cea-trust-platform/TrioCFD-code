#Script de recuperation des proprietes physiques et de la geometrie

import os, sys, math

def readPropertiesData(nomFic):
	#initialisation
	properties = {}
	properties['mu'] = -1
	properties['rho'] = -1
	properties['diffusivity'] = -1
	properties['vimp'] = -1
	properties['cimp'] = -1
	properties['diametre_jet'] = -1
	properties['velocity_decay_cte'] = -1
	properties['velocity_virtual_origin'] = -1
	properties['concentration_decay_cte'] = -1
	properties['concentration_virtual_origin'] = -1
	# ouverture des fichiers
	fic = open(nomFic,'r')
#
#
#
	for ligne in fic:
		ligne = ligne.strip().lower()
		tLigne = ligne.split()
		if ligne.find('champ_uniforme')>-1:
			if ligne.startswith('mu'):
				properties['mu'] = float(tLigne[-1])
			elif ligne.startswith('rho'):
				properties['rho'] = float(tLigne[-1])
			elif ligne.startswith('diffusivite'):
				properties['diffusivity'] = float(tLigne[-1])
		elif ligne.find('frontiere_ouverte_vitesse_imposee champ_front_uniforme')>-1:
			properties['vimp'] = float(tLigne[-1])
		elif ligne.find('frontiere_ouverte_concentration_imposee champ_front_uniforme')>-1:	
			properties['cimp'] = abs(float(tLigne[-1]))
		elif ligne.find('diametre_jet')>-1:	
			properties['diametre_jet'] = float(tLigne[-1])
		elif ligne.find('velocity_decay_cte')>-1:	
			properties['velocity_decay_cte'] = float(tLigne[-1])
		elif ligne.find('velocity_virtual_origin')>-1:	
			properties['velocity_virtual_origin'] = float(tLigne[-1])
		elif ligne.find('concentration_decay_cte')>-1:	
			properties['concentration_decay_cte'] = float(tLigne[-1])
		elif ligne.find('concentration_virtual_origin')>-1:	
			properties['concentration_virtual_origin'] = float(tLigne[-1])
#			
	fic.close()
	return properties


def ecritureFichier(properties):
	#ecriture du fichier
	nomFic = 'propertiesGeometry.dat'
	fichier = open(nomFic, 'w')

	fichier.write('%18.4f %18.2f %18.11f %18.2f %18.2f  %18.3f  %18.3f  %18.3f  %18.3f  %18.3f\n' % ( properties['mu'], properties['rho'], properties['diffusivity'], properties['vimp'], properties['cimp'], properties['diametre_jet'], properties['velocity_decay_cte'], properties['velocity_virtual_origin'], properties['concentration_decay_cte'], properties['concentration_virtual_origin'] ))
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
		properties['rho'] = float(tabLigne[ind])
		ind += 1
		properties['diffusivity'] = float(tabLigne[ind])
		ind += 1
		properties['vimp'] = float(tabLigne[ind])
		ind += 1
		properties['cimp'] = float(tabLigne[ind])
		ind += 1
		properties['diametre_jet'] = float(tabLigne[ind])
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
