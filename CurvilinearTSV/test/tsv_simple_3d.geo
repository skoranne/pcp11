SetFactory("OpenCASCADE");

// Mesh options
Mesh.CharacteristicLengthMin = 1e-07;
Mesh.CharacteristicLengthMax = 1e-05;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 10;

// Air box (margin=5e-05)
Box(1) = {-8.5e-05, -6e-05, -5e-05, 0.000232, 0.00012, 0.0001005};

Box(2) = {-3.5e-05, -1e-05, 0, 2e-05, 2e-05, 1e-07};
Box(3) = {1.5e-05, -1e-05, 0, 2e-05, 2e-05, 1e-07};
Box(4) = {-1e-05, -1e-05, 0, 2e-05, 2e-05, 1e-07};
Box(5) = {4e-05, -1e-05, 0, 2e-05, 2e-05, 1e-07};
Box(6) = {7e-05, -1e-06, 2e-07, 2e-06, 2e-06, 1e-07};
Box(7) = {9.5e-05, -1e-06, 2e-07, 2e-06, 2e-06, 1e-07};
Box(8) = {7e-05, -1e-06, 4e-07, 2e-06, 2e-06, 1e-07};
Box(9) = {9.5e-05, -1e-06, 4e-07, 2e-06, 2e-06, 1e-07};

// Fragment all volumes for conformal mesh
v() = BooleanFragments{ Volume{1, 2, 3, 4, 5, 6, 7, 8, 9}; Delete; }{};

// Identify volumes by bounding box after fragmentation
eps = 1e-12;
assigned[] = {};
layerVols_L1_0() = {};
For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  zc = (bb(2) + bb(5)) / 2;
  If (zc >= 0 - eps && zc <= 1e-07 + eps)
    layerVols_L1_0() += {v(i)};
    assigned[] += {v(i)};
  EndIf
EndFor
If (#layerVols_L1_0() > 0)
  Physical Volume(1) = {layerVols_L1_0()}; // L1_0
EndIf
layerVols_L2_0() = {};
For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  zc = (bb(2) + bb(5)) / 2;
  If (zc >= 2e-07 - eps && zc <= 3e-07 + eps)
    layerVols_L2_0() += {v(i)};
    assigned[] += {v(i)};
  EndIf
EndFor
If (#layerVols_L2_0() > 0)
  Physical Volume(2) = {layerVols_L2_0()}; // L2_0
EndIf
layerVols_L3_0() = {};
For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  zc = (bb(2) + bb(5)) / 2;
  If (zc >= 4e-07 - eps && zc <= 5e-07 + eps)
    layerVols_L3_0() += {v(i)};
    assigned[] += {v(i)};
  EndIf
EndFor
If (#layerVols_L3_0() > 0)
  Physical Volume(3) = {layerVols_L3_0()}; // L3_0
EndIf
airVols() = {};
For i In {0 : #v()-1}
  found = 0;
  For j In {0 : #assigned[]-1}
    If (assigned[j] == v(i))
      found = 1;
    EndIf
  EndFor
  If (found == 0)
    airVols() += {v(i)};
  EndIf
EndFor
If (#airVols() > 0)
  Physical Volume(4) = {airVols()}; // air
EndIf
