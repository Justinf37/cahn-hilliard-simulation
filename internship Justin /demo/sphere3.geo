// Gmsh project created on Fri May 30 10:18:40 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("boundary1", 4) = {1};
//+
Point(3) = {0, 0, 0, 0.1};
//+
Physical Point("centre", 5) = {3};
//+
Physical Curve("boundary1", 6) = {2};
