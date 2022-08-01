// Premiere maille en paroi
ypm = 0.0025;
lc = 8 * ypm;


// Diametre de la forme
R = 0.025 ;
H = 0.5 ;
l = 0.25 ;
L = 1.25 ; 
Point(1) = {-l,0,0,lc};
Point(2) = {L,0,0,lc};
Point(3) = {L,H,0,ypm};
Point(4) = {-l,H,0,ypm};

Point(5) = {-R,0,0,ypm};
Point(8) = {R,0,0,ypm};

f=1.2;
n=Ceil(Log(lc/ypm)/Log(f));
y=ypm*(1-f^(n-1))/(1-f);
Point(9) = {-l+y,H-y,0,lc};
Point(10) = {L-y,H-y,0,lc};

Line(2) = {1,5};
Line(5) = {8,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};
Line(9) = {9,10};

Point(6) = {0,0,0,ypm};
Point(7) = {0,R,0,ypm};
Circle(1) = {5,6,7};
Circle(3) = {7,6,8};
Line Loop(5) = {2,1,3,5,6,7,8,9};
Physical Line("CYLINDER") = {1,3};
Physical Line("BOTTOM") = {2,5};

Plane Surface(9) = {5};

// Frontieres
Physical Line("OUTLET") = {6};
Physical Line("TOP") = {7};
Physical Line("INLET") = {8};

// Domaine
Physical Surface("domaine") = {9};


//Field[1] = Box;
//Field[1].VIn = lc;
//Field[1].VOut = 0.001;
//Field[1].XMin = 0;
//Field[1].XMax = L;
//Field[1].YMin = 2*D;
//Field[1].YMax = H-2*D;
//Field[2] = Max;
//Field[2].FieldsList = {1};
//Background Field = 2;

Show "*";
