# Ecoulement 3D autour d'une barre de section carree 	#
# Caillaux, 04/01/99									#
# Cas test pour le maillage "Tetraedriser_Homogene"		#
#	utilisation du schema   "VEFPreP1B"					#

Mailler dom
{
	Pave Entree
	{
		Origine -5.0 -0.5 -2.0
		Longueurs 4.5 1.0 4.0
		Nombre_de_Noeuds 3 2 3
	}
	{
		Bord Entree_1 	X = -5.0 	-0.5 <= Y <= 0.5 	-2.0 <= Z <= 2.0
		Bord Entree_2 	X = -0.5 	-0.5 <= Y <= 0.5	-2.0 <= Z <= 2.0
		Bord Entree_5 	Z = -2.0	-5.0 <= X <= -0.5 	-0.5 <= Y <= 0.5
		Bord Entree_6 	Z = 2.0		-5.0 <= X <= -0.5 	-0.5 <= Y <= 0.5
	} ,

	Pave Entree_Haut
	{
		Origine -5.0 0.5 -2.0
		Longueurs 4.5 6.5 4.0
		Nombre_de_Noeuds 3 2 3
	}
	{
		Bord Entree_Haut_1 	X = -5.0	0.5 <= Y <= 7.0	 	-2.0 <= Z <= 2.0
		Bord Entree_Haut_4 	Y = 7.0		-5.0 <= X <= -0.5 	-2.0 <= Z <= 2.0
		Bord Entree_Haut_5 	Z = -2.0	-5.0 <= X <= -0.5 	0.5 <= Y <= 7.0
		Bord Entree_Haut_6 	Z = 2.0		-5.0 <= X <= -0.5 	0.5 <= Y <= 7.0
	} ,

	Pave Entree_Bas
	{
		Origine -5.0 -7.0 -2.0
		Longueurs 4.5 6.5 4.0
		Nombre_de_Noeuds 3 2 3
	}
	{
		Bord Entree_Bas_1 	X = -5.0	-7.0 <= Y <= -0.5 	-2.0 <= Z <= 2.0
		Bord Entree_Bas_3 	Y = -7.0	-5.0 <= X <= -0.5 	-2.0 <= Z <= 2.0
		Bord Entree_Bas_5 	Z = -2.0	-5.0 <= X <= -0.5 	-7.0 <= Y <= -0.5
		Bord Entree_Bas_6 	Z = 2.0		-5.0 <= X <= -0.5 	-7.0 <= Y <= -0.5
	} ,

	Pave Sortie
	{
		Origine 0.5 -0.5 -2.0
		Longueurs 15.0 1.0 4.0
		Nombre_de_Noeuds 4 2 3
	}
	{
		Bord Sortie_1 	X = 0.5 	-0.5 <= Y <= 0.5 	-2.0 <= Z <= 2.0
		Bord Sortie_2 	X = 15.5 	-0.5 <= Y <= 0.5 	-2.0 <= Z <= 2.0
		Bord Sortie_5 	Z = -2.0	-0.5 <= X <= 15.5 	-0.5 <= Y <= 0.5
		Bord Sortie_6 	Z = 2.0		-0.5 <= X <= 15.5 	-0.5 <= Y <= 0.5
	} ,

	Pave Sortie_Haut
	{
		Origine 0.5 0.5 -2.0
		Longueurs 15.0 6.5 4.0
		Nombre_de_Noeuds 4 2 3
	}
	{
		Bord Sortie_Haut_2 	X = 15.5 	0.5 <= Y <= 7.0	 	-2.0 <= Z <= 2.0
		Bord Sortie_Haut_4 	Y = 7.0		0.5 <= X <= 15.5 	-2.0 <= Z <= 2.0
		Bord Sortie_Haut_5 	Z = -2.0	-0.5 <= X <= 15.5 	0.5 <= Y <= 7.0
		Bord Sortie_Haut_6 	Z = 2.0		-0.5 <= X <= 15.5 	0.5 <= Y <= 7.0
	} ,

	Pave Sortie_Bas
	{
		Origine 0.5 -7.0 -2.0
		Longueurs 15.0 6.5 4.0
		Nombre_de_Noeuds 4 2 3
	}
	{
		Bord Sortie_Bas_2 	X = 15.5 	-7.0 <= Y <= -0.5 	-2.0 <= Z <= 2.0
		Bord Sortie_Bas_3 	Y = -7.0	0.5 <= X <= 15.5 	-2.0 <= Z <= 2.0
		Bord Sortie_Bas_5 	Z = -2.0	-0.5 <= X <= 15.5 	-7.0 <= Y <= -0.5
		Bord Sortie_Bas_6 	Z = 2.0		-0.5 <= X <= 15.5 	-7.0 <= Y <= -0.5
	} ,

	Pave Haut 
	{
		Origine -0.5 0.5 -2.0
		Longueurs 1.0 6.5 4.0
		Nombre_de_Noeuds 2 2 3
	}
	{
		Bord Haut_3 	Y = 0.5 	-0.5 <= X <= 0.5 	-2.0 <= Z <= 2.0
		Bord Haut_4 	Y = 7.0		-0.5 <= X <= 0.5 	-2.0 <= Z <= 2.0
		Bord Haut_5 	Z = -2.0	-0.5 <= X <= 0.5 	0.5 <= Y <= 7.0
		Bord Haut_6 	Z = 2.0		-0.5 <= X <= 0.5 	0.5 <= Y <= 7.0
	} ,

	Pave Bas 
	{
		Origine -0.5 -7.0 -2.0
		Longueurs 1.0 6.5 4.0
		Nombre_de_Noeuds 2 2 3
	}
	{
		Bord Bas_3 	Y = -7.0	-0.5 <= X <= 0.5 	-2.0 <= Z <= 2.0
		Bord Bas_4 	Y = -0.5 	-0.5 <= X <= 0.5 	-2.0 <= Z <= 2.0
		Bord Bas_5 	Z = -2.0	-0.5 <= X <= 0.5 	-7.0 <= Y <= -0.5
		Bord Bas_6 	Z = 2.0		-0.5 <= X <= 0.5 	-7.0 <= Y <= -0.5
	}
}
