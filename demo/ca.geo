// Gmsh project created on Tue May 27 12:22:18 2025
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
Point(101) = {0, 0, 1, 1.0};
Physical Point("fix") = {101};

// Grouper la surface pour la traction
Physical Surface("load")  = {1};

//  Grouper le volume (pas strictement nécessaire ici)
Physical Volume ("domain")= {1};

//  Forcer MSH 4.1 ASCII
Mesh.MshFileVersion = 4.1;
Mesh.Format         = 1;//+
Physical Surface(" load", 2) -= {1};
