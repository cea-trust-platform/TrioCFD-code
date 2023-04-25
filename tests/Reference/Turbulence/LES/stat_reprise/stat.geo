Mailler dom 
{
Pave Entree
	{
	   Origine 0. 1. 0.
	   Nombre_de_Noeuds 8 6 11
	   Longueurs 7. 1. 10.
	}
	{
	   Bord Entree X = 0.  1. <= Y <= 2. 0. <= Z <= 10. 
	   Bord Paroi  Y = 2.  0. <= X <= 7. 0. <= Z <= 10. 
	   Bord Paroi  Y = 1.  0. <= X <= 7. 0. <= Z <= 10. 
	   Bord Paroi  Z = 0.  0. <= X <= 7. 1. <= Y <= 2. 
	   Bord Paroi  Z = 10. 0. <= X <= 7. 1. <= Y <= 2.
	} , 

Pave Haut
	{
	   Origine 7. 1. 0.
	   Nombre_de_Noeuds 11 6 11
	   Longueurs 10. 1. 10.
	}
	{
	   Bord Paroi  Y = 2.  7. <= X <= 17. 0. <= Z <= 10. 
	   Bord Paroi  Z = 0.  7. <= X <= 17. 1. <= Y <= 2. 
	   Bord Paroi  Z = 10. 7. <= X <= 17. 1. <= Y <= 2.
	} ,

Pave SHaute
	{
	   Origine 17. 1. 0.
	   Nombre_de_Noeuds 14 6 11
	   Longueurs 13. 1. 10.
	}
	{
	   Bord SortieHaute X = 30.  1. <= Y <= 2.  0. <= Z <= 10. 
	   Bord Paroi	 Y = 2.  17. <= X <= 30. 0. <= Z <= 10. 
	   Bord Paroi	 Z = 0.  17. <= X <= 30. 1. <= Y <= 2. 
	   Bord Paroi	 Z = 10. 17. <= X <= 30. 1. <= Y <= 2.
	} ,

Pave Bas 
	{
	   Origine 7. 0. 0.
	   Nombre_de_Noeuds 11 6 11
	   Longueurs 10. 1. 10.
	}
	{
	   Bord Paroi   Y = 0.  7. <= X <= 17. 0. <= Z <= 10. 
	   Bord Paroi   X = 7.  0. <= Y <= 1.  0. <= Z <= 10. 
	   Bord Paroi   Z = 0.  7. <= X <= 17. 0. <= Y <= 1. 
	   Bord Paroi   Z = 10. 7. <= X <= 17. 0. <= Y <= 1.
	} ,

Pave SBasse
	{
	   Origine 17. 0. 0.
	   Nombre_de_Noeuds 14 6 11
	   Longueurs 13. 1. 10.
	}
	{
	   Bord SortieBasse X = 30.  0. <= Y <= 1.  0. <= Z <= 10. 
	   Bord Paroi	 Y = 0.  17. <= X <= 30. 0. <= Z <= 10. 
	   Bord Paroi	 Z = 0.  17. <= X <= 30. 0. <= Y <= 1. 
	   Bord Paroi	 Z = 10. 17. <= X <= 30. 0. <= Y <= 1.
	}
}
