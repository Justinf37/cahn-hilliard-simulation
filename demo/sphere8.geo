// Gmsh project created on Mon Jun  2 11:10:58 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0.5, 0.01};
//+
Sphere(1) = {0, 0, 0.5, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Point("centre", 4) = {1};
//+
Physical Surface("boundary1", 5) = {1};
//+
Physical Volume("volume", 6) = {1};
