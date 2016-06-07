# Maillage d une cavite 2D pour prise en compte des conditions #
#  aux limites ouvertes #
# 1- on etend la cavite hors du domaine de 2L en largeur, 1;5L au dessus #
#    et L en dessous #
# 2- La dimension L est calculee pour respecter le nombre de Rayleigh #
#    Ra= 1e+6 L=0.387315 #
# 3- on maille uniformement tout le domaine #

Mailler dom {
     Pave C1
     	  {
	  Origine 0. 0.2346
	  Nombre_de_Noeuds 40 40
	  Longueurs 0.2346 0.2346
	  }
	  {
	  bord HOT X = 0. 0.2346 <= Y <= 0.4692
	  bord ADIAB1 Y = 0.2346 0. <= X <= 0.2346
	  bord ADIAB2 Y = 0.4692 0. <= X <= 0.2346
	  } , 
     Pave C2
     	  {
	  Origine 0.2346 0.2346
	  Nombre_de_Noeuds 80 40
	  Longueurs 0.4692 0.2346
	  }
	  {
	  bord COLD3 X = 0.7038 0.2346 <= Y <= 0.4692
	  } , 
     Pave C3
     	  {
	  Origine 0.2346 0.
	  Nombre_de_Noeuds 80 40
	  Longueurs 0.4692 0.2346
	  }
	  {
	  bord COLD2 X = 0.7038  0. <= Y <= 0.2346
	  bord COLD1 Y = 0. 0.2346 <= X <= 0.7038
	  bord ADIAB3 X = 0.2346 0. <= Y <= 0.2346
	  } , 
     Pave C4
     	  {
	  Origine 0.2346 0.4692
	  Nombre_de_Noeuds 80 60
	  Longueurs 0.4692 0.3519
	  }
	  {
	  bord COLD4 X = 0.7038 0.4692 <= Y <= 0.8211
	  bord COLD5 Y = 0.8211 0.2346 <= X <= 0.7038
	  bord ADIAB4 X = 0.2346 0.4692 <= Y <= 0.8211
	  }
}

# les fichiers DOMxxxx sont crees dans le repertoire. #
# On peut mettre le chemin devant DOM si on change de repertoire. #
# Scatter DOM.Zones dom #
