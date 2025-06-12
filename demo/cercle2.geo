// Gmsh project created on Fri Jun  6 09:31:17 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
//+
Physical Surface("centre", 13) = {4};
//+
Physical Surface("boundary1", 14) = {3};
//+
Physical Volume("volume", 15) = {1};
