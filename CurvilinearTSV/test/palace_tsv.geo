SetFactory("OpenCASCADE");

// Units: microns (L0=1e-6 in Palace config)
um = 1.0;

// Geometry parameters
r = 10*um;        // TSV radius
h = 100*um;       // TSV height
pitch = 50*um;    // TSV pitch
tw = 2*um;        // trace width
tt = 2*um;        // trace thickness
tl = 30*um;       // trace length
port_gap = 2*um;  // lumped port gap
air = 80*um;      // air margin

// Air box
Box(1) = {-pitch/2-air, -air, -air, pitch+2*air, 2*air, h+2*air};

// TSV 1
Cylinder(2) = {-pitch/2, 0, 0, 0, 0, h, r};
// TSV 2
Cylinder(3) = { pitch/2, 0, 0, 0, 0, h, r};

// Top trace 1 (from TSV1 top to +x)
Box(4) = {-pitch/2, -tw/2, h, tl, tw, tt};
// Top trace 2 (from TSV2 top to +x)
Box(5) = { pitch/2, -tw/2, h, tl, tw, tt};
// Bottom trace 1
Box(6) = {-pitch/2, -tw/2, -tt, tl, tw, tt};
// Bottom trace 2
Box(7) = { pitch/2, -tw/2, -tt, tl, tw, tt};

// Port gaps (small volumes to define lumped port surfaces) - at trace ends
// Port 1: end of top trace 1
Box(8) = {-pitch/2 + tl - port_gap, -tw/2, h, port_gap, tw, tt};
// Port 2: end of top trace 2
Box(9) = { pitch/2 + tl - port_gap, -tw/2, h, port_gap, tw, tt};
// Port 3: end of bottom trace 1
Box(10) = {-pitch/2 + tl - port_gap, -tw/2, -tt, port_gap, tw, tt};
// Port 4: end of bottom trace 2
Box(11) = { pitch/2 + tl - port_gap, -tw/2, -tt, port_gap, tw, tt};

// Fragment all
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:11}; Delete; };

// Find air volume (largest)
airVol = v(0);
maxVol = 0;
For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  dx = bb(3)-bb(0); dy = bb(4)-bb(1); dz = bb(5)-bb(2);
  vol = dx*dy*dz;
  If (vol > maxVol)
    maxVol = vol;
    airVol = v(i);
  EndIf
EndFor

Physical Volume(1) = {airVol};

// Boundaries
bnd() = Boundary{ Volume{airVol}; };
gnd() = {}; pec() = {};
p1() = {}; p2() = {}; p3() = {}; p4() = {};
eps = 0.5;

For i In {0 : #bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin = bb(0); ymin = bb(1); zmin = bb(2);
  xmax = bb(3); ymax = bb(4); zmax = bb(5);
  cx = (xmin+xmax)/2; cy = (ymin+ymax)/2; cz = (zmin+zmax)/2;
  dx = xmax-xmin; dy = ymax-ymin; dz = zmax-zmin;

  onOuter = (xmin < -pitch/2-air+eps || xmax > pitch/2+air-eps ||
             ymin < -air+eps || ymax > air-eps ||
             zmin < -air+eps || zmax > h+air-eps);

  // Port surfaces: cross-sections at port gap locations (x = trace_end - port_gap)
  atP1 = (Abs(xmin - (-pitch/2 + tl - port_gap)) < eps || Abs(xmax - (-pitch/2 + tl - port_gap)) < eps) && cz > h-eps && cx < -pitch/2+tl;
  atP2 = (Abs(xmin - ( pitch/2 + tl - port_gap)) < eps || Abs(xmax - ( pitch/2 + tl - port_gap)) < eps) && cz > h-eps && cx >  pitch/2;
  atP3 = (Abs(xmin - (-pitch/2 + tl - port_gap)) < eps || Abs(xmax - (-pitch/2 + tl - port_gap)) < eps) && cz < 0+eps && cx < -pitch/2+tl;
  atP4 = (Abs(xmin - ( pitch/2 + tl - port_gap)) < eps || Abs(xmax - ( pitch/2 + tl - port_gap)) < eps) && cz < 0+eps && cx >  pitch/2;

  If (onOuter)
    gnd() += {s};
  ElseIf (atP1)
    p1() += {s};
  ElseIf (atP2)
    p2() += {s};
  ElseIf (atP3)
    p3() += {s};
  ElseIf (atP4)
    p4() += {s};
  Else
    pec() += {s};
  EndIf
EndFor

Physical Surface(2) = {gnd()};
Physical Surface(3) = {pec()};
If (#p1() > 0) Physical Surface(4) = {p1()}; EndIf
If (#p2() > 0) Physical Surface(5) = {p2()}; EndIf
If (#p3() > 0) Physical Surface(6) = {p3()}; EndIf
If (#p4() > 0) Physical Surface(7) = {p4()}; EndIf

Mesh.CharacteristicLengthMin = 2;
Mesh.CharacteristicLengthMax = 15;
Mesh.Algorithm3D = 10;
