Mailler dom 
{
Pave 1
	{
	Origine -$l 0. 0.
	Nombre_de_Noeuds $NX1 $NY1 $NZ
	Longueurs $l $H $H 
	Facteurs $FX1 $FY1 1.
	Symy
	}
	{
	Bord INLET 	X = -$l  	0. <= Y <= $H  		0. <= Z <= $H
	Bord WALL  	Y = $H   	-$l <= X <= 0. 		0. <= Z <= $H   
	Bord WALL   	Y = 0.  	-$l <= X <= 0. 		0. <= Z <= $H
	Bord SUD   	Z = 0.  	-$l <= X <= 0. 		0. <= Y <= $H  
	Bord NORD  	Z = $H  	-$l <= X <= 0. 		0. <= Y <= $H  
	} , 
Pave 2
	{
	Origine 0. 0. 0.
	Nombre_de_Noeuds $NX2 $NY2 $NZ
	Longueurs $L $H $H 
	Facteurs $FX2 $FY2 1.
	Symy	
	}
	{
	Bord OUTLET 	X = $L 		0. <= Y <= $H  		0. <= Z <= $H
	Bord WALL  	Y = $H   	0. <= X <= $L  		0. <= Z <= $H 
	Bord SUD   	Z = 0. 		0. <= X <= $L 		0. <= Y <= $H  
	Bord NORD  	Z = $H	  	0. <= X <= $L  		0. <= Y <= $H  	  
	} ,
Pave 3 
	{
	Origine 0. -$h 0.
	Nombre_de_Noeuds $NX3 $NY3 $NZ
	Longueurs $L $h $H 
	Facteurs $FX3 $FY3 1.
	Symy
	}
	{
	Bord WALL   	Y = -$h  	0. <= X <= $L  		0 <= Z <= $H 
	Bord WALL   	X = 0.  	-$h <= Y <= 0. 		0 <= Z <= $H  
	Bord OUTLET  	X = $L   	-$h <= Y <= 0. 		0 <= Z <= $H  
	Bord SUD   	Z = 0.  	0. <= X <= $L 		-$h <= Y <= 0. 
	Bord NORD  	Z = $H  	0. <= X <= $L  		-$h <= Y <= 0. 	 
	}
}
ecrire_fichier dom dom3D.geom
