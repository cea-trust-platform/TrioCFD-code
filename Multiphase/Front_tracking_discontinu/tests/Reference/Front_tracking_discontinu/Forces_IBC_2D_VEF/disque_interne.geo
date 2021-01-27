// Creation d'un cercle
Point(1) = {xc+R, yc, 0, lc};
Extrude { {0,0,1} , {xc,yc,0} , 0.5*Pi} { 
  Point{1};
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

