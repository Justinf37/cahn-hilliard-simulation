// Gmsh project created on Thu Jun  5 10:08:03 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(2) = {0, 0, 0, 0.001, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("centre", 7) = {2};
//+
Physical Surface("boundary1", 8) = {1};
//+
Physical Volume("volume", 9) = {1};
//+
Physical Surface(" centre", 7) -= {2};
//+
Physical Point("milieu", 10) = {4, 3};
