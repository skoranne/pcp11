SetFactory("OpenCASCADE");

// Units: microns (Palace L0 = 1e-6)
um = 1.0;

// Parameters
r_tsv = 5*um;          // TSV Cu radius
t_ox  = 0.5*um;        // SiO2 liner thickness
r_ox  = r_tsv + t_ox;  // outer oxide radius
h_si  = 50*um;         // substrate/TSV height
pitch = 40*um;         // GSG pitch
w_trace = 4*um;        // top trace width
t_trace = 2*um;        // top trace thickness
l_trace = 15*um;       // trace length from TSV center outward
port_gap = 2*um;       // lumped port gap
air_margin = 60*um;    // air margin
gnd_t = 2*um;          // bottom ground thickness

// Air/substrate box
Box(1) = {-pitch-air_margin, -air_margin, -gnd_t-air_margin, 2*pitch+2*air_margin, 2*air_margin, h_si+gnd_t+2*air_margin};

// Silicon substrate
Box(2) = {-pitch-air_margin, -air_margin, 0, 2*pitch+2*air_margin, 2*air_margin, h_si};

// SiO2 liners (cylinders)
Cylinder(3) = {-pitch, 0, 0, 0, 0, h_si, r_ox}; // G1
Cylinder(4) = {0,      0, 0, 0, 0, h_si, r_ox}; // S
Cylinder(5) = {pitch,  0, 0, 0, 0, h_si, r_ox}; // G2

// Cu TSVs
Cylinder(6) = {-pitch, 0, 0, 0, 0, h_si, r_tsv}; // G1
Cylinder(7) = {0,      0, 0, 0, 0, h_si, r_tsv}; // S
Cylinder(8) = {pitch,  0, 0, 0, 0, h_si, r_tsv}; // G2

// Bottom ground plane
Box(9) = {-pitch-air_margin, -air_margin, -gnd_t, 2*pitch+2*air_margin, 2*air_margin, gnd_t};

// Top ground traces connecting ground TSVs to edges (for port reference)
// Left ground trace
Box(10) = {-pitch-air_margin, -w_trace/2, h_si, air_margin + pitch - r_tsv, w_trace, t_trace};
// Right ground trace
Box(11) = {pitch + r_tsv, -w_trace/2, h_si, air_margin + pitch - r_tsv, w_trace, t_trace};
// Signal trace to left port
Box(12) = {-l_trace, -w_trace/2, h_si, l_trace - r_tsv, w_trace, t_trace};
// Signal trace to right port
Box(13) = {r_tsv, -w_trace/2, h_si, l_trace - r_tsv, w_trace, t_trace};

// Port gap volumes (between signal trace end and ground trace) - create as thin boxes
// Port 1: at left end of signal trace, gap to left ground trace
Box(14) = {-l_trace - port_gap, -w_trace/2, h_si, port_gap, w_trace, t_trace};
// Port 2: at right end of signal trace, gap to right ground trace
Box(15) = {l_trace, -w_trace/2, h_si, port_gap, w_trace, t_trace};

// Fragment all
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:15}; Delete; };

// Identify volumes by bounding box center and size heuristics
airVol() = {}; siVol() = {}; oxVol() = {}; cuVol() = {}; gndVol() = {}; traceVol() = {}; portVol() = {};

For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin = bb(0); ymin = bb(1); zmin = bb(2);
  xmax = bb(3); ymax = bb(4); zmax = bb(5);
  dx = xmax-xmin; dy = ymax-ymin; dz = zmax-zmin;
  cx = (xmin+xmax)/2; cy = (ymin+ymax)/2; cz = (zmin+zmax)/2;
  vol = dx*dy*dz;

  // Air: largest volume
  isAir = (dx > 150 && dy > 100 && dz > 100);
  isSi = (!isAir && dz > 40 && dz < 60 && dx > 150);
  isGnd = (!isAir && !isSi && zmax <= 0.5 && dz < 5);
  isOx = (!isAir && !isSi && !isGnd && Abs(cx+pitch) < 8 && dz > 40 && dx < 15);
  isOx2 = (!isAir && !isSi && !isGnd && Abs(cx) < 8 && dz > 40 && dx < 15);
  isOx3 = (!isAir && !isSi && !isGnd && Abs(cx-pitch) < 8 && dz > 40 && dx < 15);
  isCu = (!isAir && !isSi && !isGnd && (Abs(cx+pitch)<6 || Abs(cx)<6 || Abs(cx-pitch)<6) && dz > 40 && dx < 12);
  isTrace = (!isAir && !isSi && !isGnd && cz > h_si-1 && dz <= 3);
  isPort = (isTrace && (Abs(cx - (-l_trace - port_gap/2)) < 2 || Abs(cx - (l_trace + port_gap/2)) < 2));

  If (isAir)
    airVol() += {v(i)};
  ElseIf (isSi)
    siVol() += {v(i)};
  ElseIf (isGnd)
    gndVol() += {v(i)};
  ElseIf (isCu)
    cuVol() += {v(i)};
  ElseIf (isOx || isOx2 || isOx3)
    oxVol() += {v(i)};
  ElseIf (isPort)
    portVol() += {v(i)};
  ElseIf (isTrace)
    traceVol() += {v(i)};
  Else
    // Default: assign by z
    If (cz < 0)
      gndVol() += {v(i)};
    ElseIf (cz < h_si)
      siVol() += {v(i)};
    Else
      traceVol() += {v(i)};
    EndIf
  EndIf
EndFor

If (#airVol() > 0) Physical Volume(1) = {airVol()}; EndIf
If (#siVol() > 0)  Physical Volume(2) = {siVol()}; EndIf
If (#oxVol() > 0)  Physical Volume(3) = {oxVol()}; EndIf
If (#cuVol() > 0)  Physical Volume(4) = {cuVol()}; EndIf
If (#gndVol() > 0) Physical Volume(5) = {gndVol()}; EndIf
If (#traceVol() > 0) Physical Volume(6) = {traceVol()}; EndIf
// Port volumes kept in air/trace handling; we'll tag port surfaces separately

// Surface tagging
allSurfs() = Boundary{ Volume{:}; };
bndAir() = {};
For i In {0 : #airVol()-1}
  bndAir() += Boundary{ Volume{airVol(i)}; };
EndFor
// Also include surfaces adjacent to si
For i In {0 : #siVol()-1}
  bndAir() += Boundary{ Volume{siVol(i)}; };
EndFor

pecSurfs() = {}; port1Surfs() = {}; port2Surfs() = {}; absSurfs() = {};

allVols() = {};
allVols() += airVol(); allVols() += siVol(); allVols() += oxVol();
allVols() += cuVol(); allVols() += gndVol(); allVols() += traceVol(); allVols() += portVol();

bndAll() = {};
For i In {0 : #allVols()-1}
  bndAll() += Boundary{ Volume{allVols(i)}; };
EndFor

eps = 0.5;
For i In {0 : #bndAll()-1}
  s = Abs(bndAll(i));
  bb() = BoundingBox Surface{s};
  xmin = bb(0); ymin = bb(1); zmin = bb(2);
  xmax = bb(3); ymax = bb(4); zmax = bb(5);
  cx = (xmin+xmax)/2; cy = (ymin+ymax)/2; cz = (zmin+zmax)/2;
  dx = xmax-xmin; dy = ymax-ymin; dz = zmax-zmin;

  onOuter = (xmin < -pitch-air_margin+eps || xmax > pitch+air_margin-eps ||
             ymin < -air_margin+eps || ymax > air_margin-eps ||
             zmin < -gnd_t-air_margin+eps || zmax > h_si+gnd_t+air_margin-eps);

  // Port 1: x-face at -l_trace (between port gap box and signal trace)
  isP1 = (Abs(xmin - (-l_trace)) < eps || Abs(xmax - (-l_trace)) < eps) &&
         (cz > h_si-eps) && (Abs(cy) < w_trace) && (dx < eps+0.1);
  // Port 2: x-face at l_trace
  isP2 = (Abs(xmin - l_trace) < eps || Abs(xmax - l_trace) < eps) &&
         (cz > h_si-eps) && (Abs(cy) < w_trace) && (dx < eps+0.1);

  If (onOuter)
    absSurfs() += {s};
  ElseIf (isP1)
    port1Surfs() += {s};
  ElseIf (isP2)
    port2Surfs() += {s};
  Else
    pecSurfs() += {s};
  EndIf
EndFor

Physical Surface(10) = {absSurfs()};
Physical Surface(11) = {pecSurfs()};
If (#port1Surfs() > 0) Physical Surface(12) = {port1Surfs()}; EndIf
If (#port2Surfs() > 0) Physical Surface(13) = {port2Surfs()}; EndIf

Mesh.CharacteristicLengthMin = 1.0;
Mesh.CharacteristicLengthMax = 8.0;
Mesh.Algorithm3D = 10;
