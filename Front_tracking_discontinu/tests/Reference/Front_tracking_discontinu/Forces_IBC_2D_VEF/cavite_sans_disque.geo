Include "donnees.geo" ;
// Creation d'un cercle

Point(4) = {xc+2*R, yc, 0, lc};
Extrude { {0,0,1} , {xc,yc,0} , 0.5*Pi} { 
  Point{4};
}
Rotate { {0.0,0.0,1.0},{xc,yc,0}, 0.5*Pi} {
  Duplicata{ Line{1}; }
}
Rotate { {0.0,0.0,1.0},{xc,yc,0}, 0.5*Pi} {
  Duplicata{ Line{2}; }
}
Rotate { {0.0,0.0,1.0},{xc,yc,0}, 0.5*Pi} {
  Duplicata{ Line{3}; }
}
// Definition du cercle
Line Loop(1) = {1,2,3,4};


// Attention l'ordre est important: exterieur vers interieur?
Plane Surface(1) = {1};
cavite=1000;
Physical Surface(cavite)={1};

// Definition des frontieres
Physical Line("cl_1")={1,2,3,4}; 	// Cercle externe

