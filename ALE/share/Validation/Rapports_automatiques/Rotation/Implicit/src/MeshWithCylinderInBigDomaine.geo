//First cell sizes (used when points are defined):

lc = 0.0005;
lcfine = 0.0005;

//Cylinder radiuses:

Ra = 0.04;
Rb = 0.05;

//Cylinder center location:

lxc = 0.0;
lyc = lxc;

//Points definition:

Point(1) = {0, 0, 0, lcfine};
Point(2) = {Ra, 0,  0, lcfine};
Point(3) = {0, Ra, 0, lcfine};
Point(4) = {-Ra, 0, 0, lcfine};
Point(5) = {0, -Ra, 0, lcfine};
Point(6) = {Rb, 0,  0, lc};
Point(7) = {0, Rb, 0, lc};
Point(8) = {-Rb, 0, 0, lc};
Point(9) = {0, -Rb, 0, lc};

//Circles definition:
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};


//Naming boundaries:
Physical Line("CircleA") = {1,2,3,4};
Physical Line("CircleB") = {5,6,7,8};

//Lineloops:

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};

//Plane:

Plane Surface(1) = {1,2};

//Naming domain:

Physical Surface("dom") = {1} ;
