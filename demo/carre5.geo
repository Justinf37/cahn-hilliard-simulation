// Gmsh project created on Thu May 29 10:55:30 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Surface("fixed", 13) = {5};
//+
Physical Surface("Load", 14) = {2, 6, 4, 3, 1};
//+
Physical Volume("Volume", 15) = {1};
//+
Physical Surface(" Load", 14) -= {4, 1, 2, 3};
