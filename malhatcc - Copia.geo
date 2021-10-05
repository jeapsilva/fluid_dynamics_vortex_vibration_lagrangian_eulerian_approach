// Gmsh project created on Tue Apr 21 19:56:08 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {8, 0, 0, 1.0};
//+
Point(3) = {8, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Circle(1) = {4, 1, 0, 0.5, 0, 2*Pi};
//+
Line(2) = {4, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Line Loop(1) = {5, 2, 3, 4};
//+
Line Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
