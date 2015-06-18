// Premiere maille en paroi
lc = 0.1/REFINE ;
LENGTH = 0.5 ;
Point(1) = {-0.5*LENGTH,-0.5*LENGTH,0,lc};
Point(2) = {0.5*LENGTH,-0.5*LENGTH,0,lc};
Point(3) = {0.5*LENGTH,0.5*LENGTH,0,lc};
Point(4) = {-0.5*LENGTH,0.5*LENGTH,0,lc};

Line(2) = {1,2};
Line(5) = {2,3};
Line(6) = {3,4};
Line(7) = {4,1};

Line Loop(5) = {2,5,6,7};
Physical Line("BOUNDARY") = {2,5,6,7};
Plane Surface(9) = {5};
Physical Surface("domaine") = {9};

Show "*";
