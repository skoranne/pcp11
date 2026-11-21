SetFactory("OpenCASCADE");

// Dimensions (in meters)
radius = 10e-6;         // TSV radius
height = 100e-6;        // TSV height
pitch = 50e-6;          // TSV pitch (center-to-center)
traceWidth = 2e-6;      // Trace width
traceHeight = 5e-6;     // Trace thickness (from die surface)
boxSize = 300e-6;

// Air domain (large box)
Box(1) = {-boxSize/2, -boxSize/2, -boxSize/2, boxSize, boxSize, boxSize};

// TSV 1 (cylinder along Z axis) - same center as box, extruded
Cylinder(2) = {-pitch/2, 0, 0, 0, 0, height, radius};

// TSV 2 (cylinder along Z axis)
Cylinder(3) = {pitch/2, 0, 0, 0, 0, height, radius};

// Top trace 1 (rectangular prism)
Box(4) = {-pitch/2 - traceWidth/2, -traceWidth/2, traceHeight - traceHeight, traceWidth, traceWidth, traceHeight};

// Top trace 2 (rectangular prism)
Box(5) = {pitch/2 - traceWidth/2, -traceWidth/2, traceHeight - traceHeight, traceWidth, traceWidth, traceHeight};

// Bottom trace 1 (rectangular prism)
Box(6) = {-pitch/2 - traceWidth/2, -traceWidth/2, -traceHeight, traceWidth, traceWidth, traceHeight};

// Bottom trace 2 (rectangular prism)
Box(7) = {pitch/2 - traceWidth/2, -traceWidth/2, -traceHeight, traceWidth, traceWidth, traceHeight};

// Air domain (volume 1, outside all conductors)
Physical Volume(1) = {1};

// Ground - outer box surfaces (attributes 20-25 for Box 1)
Physical Surface(2) = {1, 2, 3, 4, 5, 6};

// Terminal 1: TSV 1 (cylinder surface)
Physical Surface(3) = {7, 8, 9, 10, 11, 12, 13};

// Terminal 2: TSV 2 (cylinder surface)
Physical Surface(4) = {14, 15, 16, 17, 18, 19, 20};

// Terminal 3: Top trace 1
Physical Surface(5) = {21, 22, 23, 24, 25, 26};

// Terminal 4: Top trace 2
Physical Surface(6) = {27, 28, 29, 30, 31, 32};

// Terminal 5: Bottom trace 1
Physical Surface(7) = {33, 34, 35, 36, 37, 38};

// Terminal 6: Bottom trace 2
Physical Surface(8) = {39, 40, 41, 42, 43, 44};

// Mesh size - refine near conductors
Mesh.CharacteristicLengthMin = 5e-7;
Mesh.CharacteristicLengthMax = 5e-6;
Mesh.Algorithm3D = 10; // HXT

// Transfinite mesh for better control
Mesh.Algorithm = 8;
