// Gmsh project created on Tue May 27 11:12:18 2025
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Extrude {0, 0, 1} {
  Surface{1}; 
}
//+
Physical Surface("boundary1") = {1};
//+
Physical Surface("boundary2") = {6};
//+
Physical Volume("volume") = {1};
