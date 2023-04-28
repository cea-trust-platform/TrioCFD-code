// Gmsh project created on Mon May 13 19:51:40 2013

// Pas du maillage et rayon : 
h=0.0001;
r=0.0006;

// Coordonnees du centre :
cx = 0.0;
cy = 0.0;
cz = 0.0;

// Points :
Point(1) = {cx  , cy  ,cz  , h};
Point(2) = {cx+r, cy  ,cz  , h};
Point(3) = {cx  , cy+r,cz  , h};
Point(4) = {cx  , cy  ,cz+r, h};
Point(5) = {cx-r, cy  ,cz  , h};
Point(6) = {cx  , cy-r,cz  , h};
Point(7) = {cx  , cy  ,cz-r, h};

// Arcs : 
Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,6};
Circle(4) = {6,1,2};
Circle(5) = {2,1,7};
Circle(6) = {7,1,5};
Circle(7) = {5,1,4};
Circle(8) = {4,1,2};
Circle(9) = {6,1,7};
Circle(10) = {7,1,3};
Circle(11) = {3,1,4};
Circle(12) ={4,1,6};

// Definition des contours par Line loop
Line Loop(1) = {1,11,8};
Line Loop(2) = {2,7,-11}; 
Line Loop(3) = {3,-12,-7}; 
Line Loop(4) = {4,-8,12}; 
Line Loop(5) = {5,10,-1}; 
Line Loop(6) = {-2,-10,6};
Line Loop(7) = {-3,-6,-9}; 
Line Loop(8) = {-4,9,-5}; 

// Ruled Surface : permet de définir des surfaces sphériques à partir des contours. 
Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

// Et Surface Loop : permet de définir la surface fermée engendrée par les 8 surfaces définie ci-dessous. 
Surface Loop (1) = {1,2,3,4,5,6,7,8};

// Volume final : 
Volume (1) = {1};

// Quelques commandes : 
// gmsh bulle2.geo & : Pour afficher la géométrie du résultat.
// gmsh bulle2.geo -2 : Pour mailler la surface latérale.
// gmsh bulle2.geo -3 : Pour mailler le volume de la sphère.
// gmsh bulle2.msh & : Pour afficher le maillage.
