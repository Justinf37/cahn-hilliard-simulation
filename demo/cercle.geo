// Gmsh project created on Tue May 27 13:16:10 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {-1, 0, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Point(6) = {0, 0, 1, 1.0};
//+
Point(7) = {0, 0, -1, 1.0};
//+
Circle(1) = {5, 1, 7};
//+
Circle(2) = {5, 1, 3};
//+
Circle(3) = {3, 1, 4};
//+
Circle(4) = {2, 1, 4};
//+
Circle(5) = {2, 1, 5};
//+
Circle(6) = {4, 1, 7};
//+
Circle(7) = {2, 1, 7};
//+
Circle(8) = {3, 1, 7};
//+
Circle(9) = {2, 1, 6};
//+
Circle(10) = {5, 1, 6};
//+
Circle(11) = {6, 1, 3};
//+
Circle(12) = {6, 1, 4};

//+
Line(13) = {4, 1};
//+
Line(14) = {1, 5};
//+
Line(15) = {1, 2};
//+
Line(16) = {1, 3};
//+
Line(17) = {7, 1};
//+
Line(18) = {1, 6};
//+
Curve Loop(1) = {12, -4, 9};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, -2, 10};
//+
Plane Surface(1) = {2};
//+
Curve Loop(3) = {7, -1, -5};
//+
Plane Surface(1) = {3};
