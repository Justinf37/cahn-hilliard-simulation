// Gmsh project created on Mon Jun  2 10:58:55 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, -0.5, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Point(3) = {0, 0, -0.5, 0.01};
//+
Physical Point("centre", 4) = {3};
//+
Physical Surface("boundary1", 5) = {1};
