Include "donnees.geo" ;
Include "disque_interne.geo" ;

// Creation d'un cercle
Point(4) = {xc+2*R, yc, 0, lc};
Extrude { {0,0,1} , {xc,yc,0} , 0.5*Pi} { 
  Point{4};
}
Rotate { {0.0,0.0,1.0},{xc,yc,0}, 0.5*Pi} {
  Duplicata{ Line{5}; }
}
Rotate { {0.0,0.0,1.0},{xc,yc,0}, 0.5*Pi} {
  Duplicata{ Line{6}; }
}
Rotate { {0.0,0.0,1.0},{xc,yc,0}, 0.5*Pi} {
  Duplicata{ Line{7}; }
}


// Definition du cercle
Line Loop(2) = {5,6,7,8};

// Attention l'ordre est important: exterieur vers interieur?
Plane Surface(1) = {2,1};
cavite=1000;
Physical Surface(cavite)={1};

// Definition des frontieres
Physical Line("cl_1")={1,2,3,4}; 	// Cercle interne
Physical Line("cl_2")={5,6,7,8}; 	// Cercle externe



