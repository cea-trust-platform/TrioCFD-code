Include "donnees.geo";
// Definition des 3 points de l'aile
p1 = newp;Point(p1) = {R, E, -0.5*C, lc};
p2 = newp;Point(p2) = {R, E, +0.5*C, lc} ;
p3 = newp;Point(p3) = {R+C*s, E, 0, lc} ;

// Construction de l'Aile

// 1er triangle
L1 = newl;Line(L1) = {p1,p2} ;
L2 = newl;Line(L2) = {p2,p3} ;
L3 = newl;Line(L3) = {p3,p1} ;
Line Loop(1) = {L1,L2,L3};
Plane Surface(5) = {1};
// 2eme triangle par duplication
Translate {0, -2*E, 0} { Duplicata { Surface{5}; } }
Line(4) = {p1,4};
Line(5) = {p2,5};
Line(6) = {p3,9};

// Aile
Line Loop(2) = {L2,6,-8,-5} ; Plane Surface(7) = {2};
Line Loop(3) = {L1,5,-7,-4} ; Plane Surface(8) = {3};
Line Loop(4) = {L3,4,-9,-6} ; Plane Surface(9) = {4};

// Definition du maillage surfacique (reunion des surfaces)
Aile=1000;
Physical Surface(Aile)={5,6,7,8,9};
