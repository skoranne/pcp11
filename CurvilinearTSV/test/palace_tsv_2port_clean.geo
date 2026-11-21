SetFactory("OpenCASCADE");
// Units: um, Palace L0=1e-6
// 2-port: top trace over Si with TSV to bottom ground. Port 1 top gap, Port 2 bottom gap.

tsv=10.0; tox=0.5; h=50.0; box=80.0; airTop=40.0; airBot=20.0; gndT=2.0;
traceW=4.0; traceT=2.0; traceL=25.0; gap=2.0;

Box(1) = {-box,-box,-airBot-gndT, 2*box,2*box, h+airTop+airBot+gndT}; // air
Box(2) = {-box,-box,0, 2*box,2*box,h}; // Si
Box(3) = {-tsv/2-tox,-tsv/2-tox,0, tsv+2*tox,tsv+2*tox,h}; // SiO2
Box(4) = {-tsv/2,-tsv/2,0, tsv,tsv,h}; // Cu TSV
Box(5) = {-box,-box,-gndT, 2*box,2*box,gndT}; // bottom ground
// Top metal: ground-left, port1 gap, signal trace over TSV, to right edge
Box(6) = {-box,-traceW/2,h, box-tsv/2-traceL-gap, traceW,traceT}; // left ground trace
Box(7) = {-tsv/2-traceL-gap,-traceW/2,h, gap,traceW,traceT}; // port 1 gap
Box(8) = {-tsv/2-traceL,-traceW/2,h, 2*traceL+tsv,traceW,traceT}; // signal trace spanning TSV
Box(9) = {tsv/2+traceL,-traceW/2,h, gap,traceW,traceT}; // port 2 gap (top, for symmetry) - actually make bottom port
Box(10)= {tsv/2+traceL+gap,-traceW/2,h, box-(tsv/2+traceL+gap),traceW,traceT}; // right ground trace
// Bottom: port 2 as gap between TSV bottom pad and ground
Box(11)= {-traceW/2,-traceW/2,-gndT-traceT, traceW,traceW,traceT}; // bottom pad under TSV
Box(12)= {traceW/2,-traceW/2,-gndT-traceT, traceL,traceW,traceT}; // bottom trace
Box(13)= {traceW/2+traceL,-traceW/2,-gndT-traceT, gap,traceW,traceT}; // port 2 gap
Box(14)= {traceW/2+traceL+gap,-traceW/2,-gndT-traceT, box-(traceW/2+traceL+gap),traceW,traceT}; // bottom ground ext

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:14}; Delete; };

air()={}; si()={}; ox()={}; cu()={}; gnd()={}; metal()={}; p1v()={}; p2v()={};
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  inSiZ = (zmin >= -0.1 && zmax <= h+0.1 && dz > 10);
  inBotGndZ = (zmax <= 0.1 && dz <= gndT+0.5);
  inTopZ = (zmin >= h-0.1 && zmax <= h+traceT+0.1);
  inBotTrZ = (zmin >= -gndT-traceT-0.1 && zmax <= -gndT+0.1 && dz > 0.5);
  isOuterXY = (dx > 20 || dy > 20);

  If (inBotGndZ)
    gnd() += {v(i)};
  ElseIf (inTopZ)
    isP1 = (Abs(cx - (-tsv/2-traceL-gap/2)) < 1.5 && dx < gap+1);
    isP2top = (Abs(cx - (tsv/2+traceL+gap/2)) < 1.5 && dx < gap+1);
    If (isP1) p1v() += {v(i)};
    ElseIf (isP2top) metal() += {v(i)}; // treat as metal, not a port
    Else metal() += {v(i)};
    EndIf
  ElseIf (inBotTrZ)
    isP2 = (Abs(cx - (traceW/2+traceL+gap/2)) < 1.5 && dx < gap+1);
    If (isP2) p2v() += {v(i)};
    Else metal() += {v(i)};
    EndIf
  ElseIf (inSiZ && !isOuterXY)
    If (dx <= tsv+0.2 && dy <= tsv+0.2) cu() += {v(i)};
    Else ox() += {v(i)};
    EndIf
  ElseIf (inSiZ && isOuterXY)
    si() += {v(i)};
  Else
    air() += {v(i)};
  EndIf
EndFor

If (#air() > 0) Physical Volume(1) = {air()}; EndIf
If (#si()  > 0) Physical Volume(2) = {si()}; EndIf
If (#ox()  > 0) Physical Volume(3) = {ox()}; EndIf
If (#cu()  > 0) Physical Volume(4) = {cu()}; EndIf
If (#gnd() > 0) Physical Volume(5) = {gnd()}; EndIf
If (#metal() > 0) Physical Volume(6) = {metal()}; EndIf
If (#p1v() > 0) Physical Volume(7) = {p1v()}; EndIf
If (#p2v() > 0) Physical Volume(8) = {p2v()}; EndIf

allv()={}; allv()+=air(); allv()+=si(); allv()+=ox(); allv()+=cu(); allv()+=gnd(); allv()+=metal(); allv()+=p1v(); allv()+=p2v();
bnd()={}; For i In {0:#allv()-1} bnd() += Boundary{ Volume{allv(i)}; }; EndFor
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.3;
For i In {0:#bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps || zmin < -airBot-gndT-traceT+eps || zmax > h+airTop+eps);
  isZFace = (dz < eps+0.5);
  onPortY = (Abs(cy) <= traceW/2+eps);
  onPort1X = (xmin >= -tsv/2-traceL-gap-eps && xmax <= -tsv/2-traceL+eps);
  onPort2X = (xmin >= traceW/2+traceL-eps && xmax <= traceW/2+traceL+gap+eps);
  isP1 = (isZFace && Abs(cz - h) < 0.5 && onPort1X && onPortY);
  isP2 = (isZFace && Abs(cz - (-gndT)) < 0.5 && onPort2X && onPortY);
  If (onOuter) absb() += {s};
  ElseIf (isP1) p1() += {s};
  ElseIf (isP2) p2() += {s};
  Else pec() += {s};
  EndIf
EndFor
// Fallback: search all model surfaces if ports not found
If (#p1() == 0 || #p2() == 0)
  allS() = Surface{:};
  For i In {0:#allS()-1}
    s = allS(i);
    bb() = BoundingBox Surface{s};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
    dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
    isZFace = (dz < 0.5); onPortY = (Abs(cy) <= traceW/2+eps);
    onPort1X = (xmin >= -tsv/2-traceL-gap-eps && xmax <= -tsv/2-traceL+eps);
    onPort2X = (xmin >= traceW/2+traceL-eps && xmax <= traceW/2+traceL+gap+eps);
    isP1 = (isZFace && Abs(cz - h) < 0.5 && onPort1X && onPortY);
    isP2 = (isZFace && Abs(cz - (-gndT)) < 0.5 && onPort2X && onPortY);
    If (isP1) p1() += {s}; EndIf
    If (isP2) p2() += {s}; EndIf
  EndFor
EndIf
If (#absb() > 0) Physical Surface(10) = {absb()}; EndIf
If (#pec() > 0) Physical Surface(11) = {pec()}; EndIf
If (#p1() > 0) Physical Surface(12) = {p1()}; EndIf
If (#p2() > 0) Physical Surface(13) = {p2()}; EndIf
Mesh.CharacteristicLengthMin = 2.0; Mesh.CharacteristicLengthMax = 12.0; Mesh.Algorithm3D = 10;
