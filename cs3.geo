SetFactory("OpenCASCADE");

// Dimensions
w = 3; h = 5; d = 2; // Width, Height, Depth
gap = 5;

// Create 3 Rectangular Prisms (Conductors)
Box(1) = {0, 0, 0, w, h, d};
Box(2) = {w+gap, 0, 0, w, h, d};
Box(3) = {2*(w+gap), 0, 0, w, h, d};

// Create a large surrounding box (Dielectric/Air)
Box(4) = {-10, -10, -10, 50, 30, 20};

// Fragment to ensure matching meshes at interfaces
BooleanFragments{ Volume{4}; Delete; }{ Volume{1,2,3}; Delete; }

// Define Physical Groups (3D Volumes)
Physical Volume("Rect1", 10) = {1};
Physical Volume("Rect2", 11) = {2};
Physical Volume("Rect3", 12) = {3};
Physical Volume("Diel", 13) = {4};
