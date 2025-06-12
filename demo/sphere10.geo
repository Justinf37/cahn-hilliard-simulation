// ------------------------------------------------------------
// fichier : sphere_center_corrected.geo
// ------------------------------------------------------------
SetFactory("OpenCASCADE");
// 1 Paramètres de base
lc = 0.0000001;       // taille de maille souhaitée (à ajuster)
R  = 1.0;       // rayon de la sphère

// 2 Définition explicite d'un POINT au centre (0,0,0)
//    Ici l'ID du point est "1", et on lui donne comme coord. (0,0,0).
//    λc (= lc) indique la « taille locale » de maillage autour de ce point,
//    mais ce n'est pas critique pour qu'il soit créé : seul l'Embed compte.
Point(1) = { 0.0, 0.0, 0.0, lc };

// 3 Création du VOLUME sphérique, CENTRÉ à (0,0,0), de rayon R = 1
//    On utilise la primitive Sphere (Gmsh 4.x). Son ID est "2".
Sphere(2) = { 0.0, 0.0, 0.0, R, lc };

// 4 On FORCE l’inclusion du POINT #1 (à l’origine) dans le VOLUME #2
//    Sans cette ligne Mesh.Embed, Gmsh ne mettra pas obligatoirement
//    un nœud au centre, et vous risqueriez de voir un autre point (“plus ou
//    moins au hasard”) comme nœud intérieur.
//    Ici, « 0D » signifie que l’on embedde une entité de dimension 0 (un point),
//    et « Volume {2} » indique le cube de référence (id = 2) où l’on veut
//    que le nœud apparaisse.
Mesh.Embed 0D { 1 } In Volume { 2 };

// 5 (Optionnel mais très conseillé) Étiqueter ce point pour le repérer facilement
//    dans le .msh avec le nom “Center”. Lorsque vous ouvrirez le .msh dans
//    ParaView ou dans GridapGmsh, vous aurez un Physical Group 0D “Center”.
//    Cela permet de récupérer directement l’indice du nœud sans tolérance.
Physical Point("centre") = { 1 };

// 6 (Optionnel) Si vous voulez taguer la surface extérieure pour poser
//    une condition de Dirichlet u=0 sur toute la coque, par exemple :
//    Ci-dessous, on regroupe TOUTES les surfaces générées pour la sphère
//    dans un Physical Group nommé “Boundary”.
//    Gmsh détecte automatiquement quelles surfaces appartiennent au volume 2.
Physical Surface("boundary1") = { any };

// 7 (Optionnel) Étiqueter le volume lui-même (le « domaine ») sous un tag :
//    Utile si vous voulez, plus tard, récupérer le domaine complet dans Gridap.
Physical Volume("Omega") = { 2 };
