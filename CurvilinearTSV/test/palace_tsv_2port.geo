SetFactory("OpenCASCADE");
// Units: um, Palace L0=1e-6
// 2-port TSV via for validation: Si substrate, SiO2 liner, Cu TSV, PEC grounds.
// Port 1: top, Port 2: bottom. Lumped ports across gaps.

tsv = 10.0;    // Cu TSV side
tox = 0.5;     // SiO2 liner
h   = 50.0;    // Si height
box = 80.0;    // half lateral extent
airTop = 40.0;
airBot = 40.0;
gndT = 2.0;
traceW = 4.0;
traceT = 2.0;
traceL = 30.0;
gap = 2.0;

// 1 Air
Box(1) = {-box,-box,-airBot-gndT, 2*box,2*box, h+airTop+airBot+gndT};
// 2 Si
Box(2) = {-box,-box,0, 2*box,2*box,h};
// 3 SiO2
Box(3) = {-tsv/2-tox,-tsv/2-tox,0, tsv+2*tox,tsv+2*tox,h};
// 4 Cu TSV
Box(4) = {-tsv/2,-tsv/2,0, tsv,tsv,h};
// 5 Bottom ground plane
Box(5) = {-box,-box,-gndT, 2*box,2*box,gndT};
// 6 Top ground plane
Box(6) = {-box,-box,h, 2*box,2*box,gndT};

// Top signal path: pad on TSV, trace going +x, port gap, ground extension to edge
// 7 Top signal pad
Box(7) = {-traceW/2,-traceW/2,h, traceW,traceW,traceT};
// 8 Top signal trace
Box(8) = {traceW/2,-traceW/2,h, traceL-traceW/2,traceW,traceT};
// 9 Top port gap
Box(9) = {traceL,-traceW/2,h, gap,traceW,traceT};
// 10 Top ground extension after gap
Box(10)= {traceL+gap,-traceW/2,h, box-(traceL+gap),traceW,traceT};

// Bottom signal path: pad under TSV, trace going +x, port gap, ground extension
// 11 Bottom signal pad
Box(11)= {-traceW/2,-traceW/2,-gndT-traceT, traceW,traceW,traceT};
// 12 Bottom signal trace
Box(12)= {traceW/2,-traceW/2,-gndT-traceT, traceL-traceW/2,traceW,traceT};
// 13 Bottom port gap
Box(13)= {traceL,-traceW/2,-gndT-traceT, gap,traceW,traceT};
// 14 Bottom ground extension after gap
Box(14)= {traceL+gap,-traceW/2,-gndT-traceT, box-(traceL+gap),traceW,traceT};

// Connect TSV to top/bottom pads with small via pads inside ground planes
// 15 Top connection through top ground plane opening (signal pad area cutout handled by fragment)
// 16 Bottom connection through bottom ground plane opening

// Fragment
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:14}; Delete; };

// Classify
air()={}; si()={}; ox()={}; cu()={}; pecGnd()={}; pecMetal()={}; port1v()={}; port2v()={};
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  inSiZ = (zmin >= -0.1 && zmax <= h+0.1 && dz > 10);
  inTopGndZ = (zmin >= h-0.1 && zmax <= h+gndT+0.1 && dz > 0.5 && dz <= gndT+0.5);
  inBotGndZ = (zmin >= -gndT-0.1 && zmax <= 0.1 && dz > 0.5 && dz <= gndT+0.5);
  inTopTraceZ = (zmin >= h-0.1 && zmax <= h+traceT+0.1 && dz > 0.5);
  inBotTraceZ = (zmin >= -gndT-traceT-0.1 && zmax <= -gndT+0.1 && dz > 0.5);
  isOuterXY = (dx > 20 || dy > 20);

  If (inSiZ && !isOuterXY)
    If (dx <= tsv+0.2 && dy <= tsv+0.2) cu() += {v(i)};
    Else ox() += {v(i)};
    EndIf
  ElseIf (inSiZ && isOuterXY)
    si() += {v(i)};
  ElseIf (inTopTraceZ && !isOuterXY && Abs(cx-traceL-gap/2) < 1.5 && Abs(cy) < traceW/2+0.5)
    port1v() += {v(i)};
  ElseIf (inBotTraceZ && !isOuterXY && Abs(cx-traceL-gap/2) < 1.5 && Abs(cy) < traceW/2+0.5)
    port2v() += {v(i)};
  ElseIf (inTopGndZ || inBotGndZ)
    pecGnd() += {v(i)};
  ElseIf (inTopTraceZ || inBotTraceZ)
    pecMetal() += {v(i)};
  Else
    air() += {v(i)};
  EndIf
EndFor

If (#air() > 0) Physical Volume(1) = {air()}; EndIf
If (#si()  > 0) Physical Volume(2) = {si()}; EndIf
If (#ox()  > 0) Physical Volume(3) = {ox()}; EndIf
If (#cu()  > 0) Physical Volume(4) = {cu()}; EndIf
If (#pecGnd() > 0) Physical Volume(5) = {pecGnd()}; EndIf
If (#pecMetal() > 0) Physical Volume(6) = {pecMetal()}; EndIf
If (#port1v() > 0) Physical Volume(7) = {port1v()}; EndIf
If (#port2v() > 0) Physical Volume(8) = {port2v()}; EndIf

allv()={}; allv()+=air(); allv()+=si(); allv()+=ox(); allv()+=cu(); allv()+=pecGnd(); allv()+=pecMetal(); allv()+=port1v(); allv()+=port2v();
bnd()={}; For i In {0:#allv()-1} bnd() += Boundary{ Volume{allv(i)}; }; EndFor

pec()={}; p1()={}; p2()={}; absb()={};
eps=0.3;
For i In {0:#bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps || zmin < -airBot-gndT-traceT+eps || zmax > h+airTop+eps);
  onPort1Z = (zmin >= h-eps && zmax <= h+traceT+eps);
  onPort2Z = (zmin >= -gndT-traceT-eps && zmax <= -gndT+eps);
  onPortY = (Abs(cy) <= traceW/2+eps);
  isXFace = (dx < eps+0.1);
  isP1 = (isXFace && Abs(cx-traceL) < eps && onPort1Z && onPortY);
  isP2 = (isXFace && Abs(cx-traceL) < eps && onPort2Z && onPortY);
  If (onOuter) absb() += {s};
  ElseIf (isP1) p1() += {s};
  ElseIf (isP2) p2() += {s};
  Else pec() += {s};
  EndIf
EndFor

If (#absb() > 0) Physical Surface(10) = {absb()}; EndIf
If (#pec() > 0) Physical Surface(11) = {pec()}; EndIf
If (#p1() > 0) Physical Surface(12) = {p1()}; EndIf
If (#p2() > 0) Physical Surface(13) = {p2()}; EndIf

Mesh.CharacteristicLengthMin = 2.0;
Mesh.CharacteristicLengthMax = 12.0;
Mesh.Algorithm3D = 10;
