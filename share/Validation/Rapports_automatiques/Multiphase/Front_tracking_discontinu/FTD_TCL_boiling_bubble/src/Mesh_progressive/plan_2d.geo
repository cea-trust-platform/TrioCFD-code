// L : size of domaine in [m]
L = 1. ;
// N : Number of nodes if No refined are used 
N = 20; 
// FN : refined factor in the refined zone
FN = 10; 
// dx (dy) : Percentage of refined domaine in x (y)-dir [-]
dx = 0.2; 
dy = 0.2; 




// refined zone in the bottom
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L*(1-dx)/2., 0, 0, 1.0};
Point(3) = {L*(1+dx)/2., 0, 0, 1.0};
Point(4) = {L, 0, 0, 1.0};
Point(5) = {L, L*dy, 0, 1.0};
Point(6) = {L, L, 0, 1.0};
Point(7) = {L*(1+dx)/2., L, 0, 1.0};
Point(8) = {L*(1-dx)/2., L, 0, 1.0};
Point(9) = {0, L, 0, 1.0};
Point(10) = {0, L*dy, 0, 1.0};
Point(11) = {L*(1-dx)/2., L*dy, 0, 1.0};
Point(12) = {L*(1+dx)/2., L*dy, 0, 1.0};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};

Line(11) = {5, 12};
Line(12) = {12, 11};
Line(13) = {11, 10};
Line(14) = {11, 8};
Line(15) = {12, 7};
Line(16) = {3, 12};
Line(17) = {2, 11};

Line Loop(18) = {1, 17, 13, 10}; 
Line Loop(19) = {-12, 15, 7, -14}; 
Line Loop(23) = {2, 16, 12, -17}; 
Line Loop(20) = {3, 4, 11, -16}; 
Line Loop(21) = {-11, 5, 6, -15}; 
Line Loop(22) = {-13, 14, 8, 9};	

Plane Surface (18) = {18}; 
Plane Surface (19) = {19};
Plane Surface (20) = {20};
Plane Surface (21) = {21};
Plane Surface (22) = {22};
Plane Surface (23) = {23};

// creat mesh point
Transfinite Line {-1, 3, 13, -11, 8, -6} = Ceil(N*(1.-dx)/2.) Using Progression 1.2;
Transfinite Line {4, 16, 17, 10} = Ceil(N*dy*FN) Using Progression 1;
Transfinite Line {2, 12, 7} = Ceil(N*dx*FN) Using Progression 1;
Transfinite Line {5, -9, 14, 15} = Ceil(N*(1.-dy)*FN) Using Progression 1;


Transfinite Surface {18, 19, 20, 21, 22, 23} ;
Recombine Surface {18, 19, 20, 21, 22, 23};

Physical Line("left") = {9, 10};
Physical Line("top") = {6, 7, 8};
Physical Line("right") = {4, 5};
Physical Line("bottom") = {1, 2, 3};


