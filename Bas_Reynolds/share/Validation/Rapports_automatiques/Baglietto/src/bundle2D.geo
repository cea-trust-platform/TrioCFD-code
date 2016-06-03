diameter = 0.12;
pitch = 0.1404;

R = diameter/2;
P = (pitch-diameter)/2.;

xc = 0.;
yc = 0.;

theta = 2*Pi/6; 

Nx = 10;
Ny = 10;

x0 = xc; y0 = yc;
x1 = xc; y1 = yc+R;
x2 = xc; y2 = yc+R+P;
x3 = xc+(R+P)/Tan(theta); y3 = yc+R+P;
x4 = xc+R*Cos(theta);     y4 = yc+R*Sin(theta);

dist = Sqrt((x4-x3)*(x4-x3)+(y4-y3)*(y4-y3));

Printf("(x1, y1) = (%g, %g)",x1,y1);
Printf("(x2, y2) = (%g, %g)",x2,y2);
Printf("(x3, y3) = (%g, %g)",x3,y3);
Printf("(x4, y4) = (%g, %g)",x4,y4);


dx = Sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2))/(Nx-1);

Point(0) = {x0,y0,0.,dx};
Point(1) = {x1,y1,0.,dx};
Point(2) = {x2,y2,0.,dx};
Point(3) = {x3,y3,0.,dx};
Point(4) = {x4,y4,0.,dx};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Circle(4) = {4,0,1};

Line Loop(5) = {1,2,3,4};

Plane Surface(6) = {5};

//Transfinite Line{-1,3} = Nx Using Progression 1.1;
//Transfinite Line{1,-3} = Nx Using Progression 1.04;
Transfinite Line{1,3} = Nx;
Transfinite Line{2,4} = Ny;

Transfinite Surface{6} = {1, 2, 3, 4};

//Mesh.Smoothing = 100;

Printf("4*dx = %g",4*dx);

Physical Line("paroi_rod") = {4};
Physical Line("symhaut")   = {2};
Physical Line("symdroite") = {3};
Physical Line("symgauche") = {1};
Physical Surface("Face1")  = {6};
