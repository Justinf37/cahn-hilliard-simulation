// Gmsh project created on Tue Jun  3 10:11:59 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
//+
//+
Point(1000) = {0, 0, 0, 0.00001};  // Un point au centre
Physical Point("centre", 1001) = {1000};
Point(1001) = {0.0001, 0, 0, 0.00001};  // proche, mais différent
Line(999) = {1000, 1001};  // sert à ce que Gmsh conserve le point
Physical Volume("volume", 13) = {1};
//+
Physical Surface("boundary1", 14) = {2, 1, 6, 4, 5, 3};

