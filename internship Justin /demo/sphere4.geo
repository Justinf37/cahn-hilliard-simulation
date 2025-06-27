// Gmsh project created on Fri Jun  6 09:49:46 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1.0, -Pi/2, Pi/2, 2*Pi};
//+
Point(3) = {0, 0, 0, 0.000001};
//+
Physical Point("centre", 4) = {3};
//+
Physical Surface("boundary1", 5) = {1};
//+
Physical Volume("volume", 6) = {1};
//+
Physical Curve("boundary1", 7) = {2};
