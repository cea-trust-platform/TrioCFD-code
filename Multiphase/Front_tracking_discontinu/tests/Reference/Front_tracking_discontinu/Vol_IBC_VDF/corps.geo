Include "donnees.geo";
// Creation d'un cercle
Point(1) = {R, 0, -0.5*L, lc};
Extrude { {0,0,1} , {0,0,0} , 0.5*Pi} { 
  Point{1};
}
Rotate { {0.0,0.0,1.0},{0.0,0.0,0.0}, 0.5*Pi} {
  Duplicata{ Line{1}; }
}
Rotate { {0.0,0.0,1.0},{0.0,0.0,0.0}, 0.5*Pi} {
  Duplicata{ Line{2}; }
}
Rotate { {0.0,0.0,1.0},{0.0,0.0,0.0}, 0.5*Pi} {
  Duplicata{ Line{3}; }
}
// Surface du cercle
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Extrusion pour avoir le cylindre
Extrude {0, 0, L} { Surface{1}; }

// Definition du maillage surfacique (reunion des surfaces)
Corps=1000;
Physical Surface(Corps)={1,13,17,21,25,26};
