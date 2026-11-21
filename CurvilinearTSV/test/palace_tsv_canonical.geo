SetFactory("OpenCASCADE");

// Units: microns (Palace L0=1e-6)
// Canonical 2-port single TSV through silicon with SiO2 liner and bottom ground.
// Box-based for robust OpenCASCADE meshing.

tsv_w = 10.0;    // Cu TSV square side
t_ox  = 0.5;     // SiO2 liner thickness
h_si  = 50.0;    // Si thickness / TSV height
box   = 50.0;    // lateral half-extent from TSV center
air_top = 30.0;
air_bot = 10.0;
gnd_t = 2.0;     // bottom ground plane thickness
trace_w = 4.0;
trace_t = 2.0;
trace_l = 20.0;  // trace length outward
port_gap = 2.0;

// Domain
Box(1) = {-box, -box, -air_bot-gnd_t, 2*box, 2*box, h_si + air_top + air_bot + gnd_t};
// Silicon
Box(2) = {-box, -box, 0, 2*box, 2*box, h_si};
// SiO2 liner
Box(3) = {-tsv_w/2-t_ox, -tsv_w/2-t_ox, 0, tsv_w+2*t_ox, tsv_w+2*t_ox, h_si};
// Cu TSV
Box(4) = {-tsv_w/2, -tsv_w/2, 0, tsv_w, tsv_w, h_si};
// Bottom PEC ground plane
Box(5) = {-box, -box, -gnd_t, 2*box, 2*box, gnd_t};
// Top signal trace (+x)
Box(6) = {tsv_w/2, -trace_w/2, h_si, trace_l, trace_w, trace_t};
// Top ground trace (-x) connects to outer boundary reference
Box(7) = {-box, -trace_w/2, h_si, box - tsv_w/2, trace_w, trace_t};
// Port gap between signal trace end and ground extension
Box(8) = {tsv_w/2 + trace_l, -trace_w/2, h_si, port_gap, trace_w, trace_t};
// Ground extension after port
Box(9) = {tsv_w/2 + trace_l + port_gap, -trace_w/2, h_si, box - (tsv_w/2 + trace_l + port_gap), trace_w, trace_t};

// Fragment
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:9}; Delete; };

// Classify volumes by bounding box
air() = {}; si() = {}; ox() = {}; cu() = {}; gnd() = {}; metal() = {}; portvol() = {};
For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;

  inSiZ = (zmin >= -0.1 && zmax <= h_si + 0.1 && dz > 10);
  inGndZ = (zmax <= 0.1);
  inTopMetalZ = (zmin >= h_si - 0.1 && zmax <= h_si + trace_t + 0.1);
  isOuterXY = (dx > 20 || dy > 20);

  If (inGndZ && dz <= gnd_t + 0.1)
    gnd() += {v(i)};
  ElseIf (inTopMetalZ)
    isPort = (Abs(cx - (tsv_w/2 + trace_l + port_gap/2)) < 1.5 && Abs(cy) < trace_w+0.1 && dx < port_gap+1.0);
    If (isPort)
      portvol() += {v(i)};
    Else
      metal() += {v(i)};
    EndIf
  ElseIf (inSiZ && !isOuterXY)
    // Central TSV stack: Cu core or oxide
    If (dx <= tsv_w + 0.2 && dy <= tsv_w + 0.2)
      cu() += {v(i)};
    Else
      ox() += {v(i)};
    EndIf
  ElseIf (inSiZ && isOuterXY)
    si() += {v(i)};
  Else
    air() += {v(i)};
  EndIf
EndFor

If (#air() > 0) Physical Volume(1) = {air()}; EndIf
If (#si()  > 0) Physical Volume(2) = {si()};  EndIf
If (#ox()  > 0) Physical Volume(3) = {ox()};  EndIf
If (#cu()  > 0) Physical Volume(4) = {cu()};  EndIf
If (#gnd() > 0) Physical Volume(5) = {gnd()}; EndIf
If (#metal() > 0) Physical Volume(6) = {metal()}; EndIf
If (#portvol() > 0) Physical Volume(7) = {portvol()}; EndIf

// Surfaces
allv() = {};
allv() += air(); allv() += si(); allv() += ox(); allv() += cu();
allv() += gnd(); allv() += metal(); allv() += portvol();
bnd() = {};
For i In {0 : #allv()-1}
  bnd() += Boundary{ Volume{allv(i)}; };
EndFor

pec() = {}; port1() = {}; port2() = {}; absorbing() = {};
eps = 0.3;
For i In {0 : #bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;

  onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps ||
             zmin < -air_bot-gnd_t+eps || zmax > h_si+air_top+eps);

  // Port surfaces: x-faces of port gap volume
  xPortLeft = tsv_w/2 + trace_l;
  xPortRight = tsv_w/2 + trace_l + port_gap;
  onPortZ = (zmin >= h_si-eps && zmax <= h_si+trace_t+eps);
  onPortY = (Abs(cy) <= trace_w/2 + eps);
  isXFace = (dx < eps+0.1);
  isPort1 = (isXFace && Abs(cx-xPortLeft) < eps && onPortZ && onPortY);
  isPort2 = (isXFace && Abs(cx-xPortRight) < eps && onPortZ && onPortY);

  If (onOuter)
    absorbing() += {s};
  ElseIf (isPort1)
    port1() += {s};
  ElseIf (isPort2)
    port2() += {s};
  Else
    pec() += {s};
  EndIf
EndFor

If (#absorbing() > 0) Physical Surface(10) = {absorbing()}; EndIf
If (#pec() > 0) Physical Surface(11) = {pec()}; EndIf
If (#port1() > 0) Physical Surface(12) = {port1()}; EndIf
If (#port2() > 0) Physical Surface(13) = {port2()}; EndIf

Mesh.CharacteristicLengthMin = 1.5;
Mesh.CharacteristicLengthMax = 8.0;
Mesh.Algorithm3D = 10;
