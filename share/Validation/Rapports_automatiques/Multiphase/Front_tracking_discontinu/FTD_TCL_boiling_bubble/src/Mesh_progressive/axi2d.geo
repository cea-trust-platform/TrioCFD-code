// L : size of domaine in [m]
LX = 0.003 ;
LY = 0.006 ;
// N : Number of nodes if No refined are used 
NX = 20; 
NY = 40;
// FN : refined factor in the refined zone
FN = 10; 
// dr : Percentage of refined domaine in r-dir [-]
dr = 0.1; 


// refined zone in the bottom
Point(1) = {0, 0, 0, 1.0};
Point(2) = {LX*dr, 0, 0, 1.0};
Point(3) = {LX, 0, 0, 1.0};
Point(4) = {LX, LY, 0, 1.0};
Point(5) = {LX*dr, LY, 0, 1.0};
Point(6) = {0, LY, 0, 1.0};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {5, 2};



// creat mesh point
Transfinite Line {1, 5} = Ceil(NX*dr*FN) Using Progression 1;
Transfinite Line {3, 6, 7} = Ceil(NY*FN) Using Progression 1;
Transfinite Line {2, -4} = Ceil(NX*(1.-dr)) Using Progression 1.2;

Line Loop(8) = {1, -7, 5, 6}; 
Line Loop(9) = {2, 3, 4, 7}; 
Plane Surface (8) = {8}; 
Plane Surface (9) = {9}; 

Transfinite Surface {8, 9} ;
Recombine Surface {8, 9};

Physical Line("left") = {6};
Physical Line("top") = {4, 5};
Physical Line("right") = {3};
Physical Line("bottom") = {1, 2};
// Physical Surface("dom") = {9};


