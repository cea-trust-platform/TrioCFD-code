//Domain dimensions:

lx = 10.0; 
ly = 5.0;

//First cell size (used when points are defined):

lc = 0.05;

//Square side length:

D = 1.0;

//Square center location:

lxc = 1.5;
lyc = 2.5;

//Points definition:

Point(1) = {0, 0, 0, lc};
Point(2) = {lx, 0,  0, lc};
Point(3) = {lx, ly, 0, lc};
Point(4) = {0, ly, 0, lc};
Point(5) = {lxc-0.5*D, lyc-0.5*D, 0, lc};
Point(6) = {lxc+0.5*D, lyc-0.5*D, 0, lc};
Point(7) = {lxc+0.5*D, lyc+0.5*D, 0, lc};
Point(8) = {lxc-0.5*D, lyc+0.5*D, 0, lc};

//Lines definition:
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

//Naming boundaries:
Physical Line("Wall") = {1,2,3,4};
Physical Line("SquareSouth") = {5};
Physical Line("SquareEast") = {6};
Physical Line("SquareNorth") = {7};
Physical Line("SquareWest") = {8};

//Lineloops:

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};

//Plane:

Plane Surface(1) = {1,2};

//Naming domain:

Physical Surface("dom") = {1} ;
