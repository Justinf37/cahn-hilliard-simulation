SetFactory("OpenCASCADE");

r = 1.0;
Sphere(1) = {0,0,0,r};

// Crée une surface entière et un point central
Point(100) = {0,0,0, 0.1};

Physical Surface("boundary1") = {1};
Physical Point("centre") = {100};
