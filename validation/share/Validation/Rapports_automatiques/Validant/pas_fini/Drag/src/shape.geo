// Variables definition
lc = 0.02;
// First cell size (used when points are defined)
lc1 = lc * 8;
// Second cell size
lc2 = lc / 2;

// Circle diameter
D = 0.14 ;
E = D ;
param = 1;
H = param * 10 * D ;
X = param * 5 * D ;
L = param * 10 * D + X + E;

// Points definition 
Point(1) = {0,0,0,lc1};
Point(2) = {L,0,0,lc1};
Point(3) = {L,H,0,lc1};
Point(4) = {0,H,0,lc1};
Point(5) = {X,0,0,lc2};
Point(8) = {X+E,0,0,lc2};

// Lines definition
Line(2) = {1,5};
Line(5) = {8,2};
Line(6) = {3,2};
Line(7) = {3,4};
Line(8) = {4,1};


// 1/4 Circle definition
Point(6) = {X+D/2,0,0,lc2}; // Center
Point(7) = {X+D/2,D/2,0,lc2};
// 3 points for the circle arc
Circle(1) = {5,6,7};
Line(3) = {7,6};
Line(4) = {6,8};
// A circle arc is STRICTLY smaller than Pi


// Naming the boundaries
Physical Line("Shape") = {1,3};
Physical Line("Axis") = {2,4,5};
Physical Line("Outlet") = {6};
Physical Line("Top") = {7};
Physical Line("Inlet") = {8};

// A lineloop is a loop on several lines
// for defining/orienting a surface
// Use negative lines to reverse the
// orientation of the line
Line Loop(1) = {2,1,3,4,5,-6,7,8};
/// The surfac will use the lineloop
Plane Surface(1) = {1};

// Naming the domain
Physical Surface("domain") = {1};
