// Use OpenCASCADE for easier primitive definition
SetFactory("OpenCASCADE");

// Define variables for easy adjustment
w = 3;   // Width
h = 5;   // Height
s = 5;   // Spacing (Distance between start of each rectangle)
lc = 0.5; // Characteristic length (mesh size)

// Create the three rectangles
// Rectangle(ID) = {x_min, y_min, z_min, width, height};
Rectangle(1) = {0, 0, 0, w, h};
Rectangle(2) = {s, 0, 0, w, h};
Rectangle(3) = {5*s, 0, 0, w, h};
Rectangle(4) = {-5, -5, 0, 50, 15};
//out[] = BooleanFragments{ Surface{box}; Delete; }{ Surface{r1, r2, r3}; Delete; };
// 3. Boolean Fragments
// This "cuts" the rectangles into the box so they share boundaries.
// 'Delete' ensures we don't have overlapping duplicate surfaces.
BooleanFragments{ Surface{4}; Delete; }{ Surface{1,2,3}; Delete; }

// 4. Physical Regions
// Note: After Boolean operations, IDs can change. 
// Using the original IDs usually works if 'Delete' is handled correctly.
Physical Surface("Rect1", 10) = {1};
Physical Surface("Rect2", 11) = {2};
Physical Surface("Rect3", 12) = {3};
Physical Surface("Dielectric", 13) = {4};

// We need to define the boundaries for capacitance (Q = CV and Q is on the boundary)
Physical Line("BndRect1", 110) = Boundary{ Surface{1}; };
Physical Line("BndRect2", 111) = Boundary{ Surface{2}; };
Physical Line("BndRect3", 113) = Boundary{ Surface{3}; };

// Define the mesh size for all points
Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;

// Optional: Synchronize to ensure the CAD engine 
// communicates with the Gmsh mesh generator
Coherence;
