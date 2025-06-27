// Gmsh project created on Thu May 29 11:19:03 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("Load", 4) = {1};

//+
Point(3) = {0, -0, 0, 0.1};
//+
Physical Point("fixed", 6) = {3};
//+
Physical Volume("volume", 7) = {1};
