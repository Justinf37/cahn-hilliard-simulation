// sphere.geo
// ------------------------------
// 1. On utilise le kernel OpenCASCADE pour la primitive "Sphere"
SetFactory("OpenCASCADE");

// 2. Définition de la sphère
r = 0.5;
Sphere(1) = {0.0, 0.0, 0.0, r};

// 3. Physical group pour la surface extérieure
//    La primitive Sphere(1) génère automatiquement une Surface de tag 1
Physical Surface("boundary1") = {1};

// 4. Point physique au centre (0,0,0)
//    On définit d’abord un point avec un maillage très fin local
Point(2) = {0.0, 0.0, 0.0, r/10};
//    Puis on marque ce point
Physical Point("centre") = {2};

// 5. Maillage
Mesh.CharacteristicLengthMax = 0.1; // ajustez pour raffiner/relâcher
Mesh 3;                             // génération 3D