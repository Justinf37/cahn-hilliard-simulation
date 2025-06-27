// Gmsh project created on Fri Jun  6 09:48:41 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.1, -Pi/2, Pi/2, 2*Pi};
//+
Point(3) = {-0, 0, -0, 0.1};
