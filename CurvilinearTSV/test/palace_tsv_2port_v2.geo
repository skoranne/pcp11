SetFactory("OpenCASCADE");
// Units: um, Palace L0=1e-6
// 2-port TSV: Cu TSV through SiO2/Si, PEC ground planes with signal cutouts.
// Port 1: top gap (+X), Port 2: bottom gap (+X). Ports on bottom face of gap.

tsv=10.0; tox=0.5; h=50.0; box=80.0; airTop=40.0; airBot=20.0; gndT=2.0;
traceW=6.0; traceT=2.0; traceL=25.0; gap=3.0; padW=12.0; // pad around TSV

Box(1) = {-box,-box,-airBot-gndT, 2*box,2*box, h+airTop+airBot+gndT}; // air
Box(2) = {-box,-box,0, 2*box,2*box,h}; // Si
Box(3) = {-tsv/2-tox,-tsv/2-tox,0, tsv+2*tox,tsv+2*tox,h}; // SiO2
Box(4) = {-tsv/2,-tsv/2,0, tsv,tsv,h}; // Cu TSV

// Ground planes with holes cut for signal paths
Box(5) = {-box,-box,-gndT, 2*box,2*box,gndT}; // bottom ground full
Box(6) = {-box,-box,h, 2*box,2*box,gndT};      // top ground full

// Top signal path cutout in top ground: from pad around TSV to right edge
Box(7) = {-padW/2,-padW/2,h, padW/2 + traceL+gap, padW, gndT}; // cutout strip to port gap
BooleanDifference{ Volume{6}; }{ Volume{7}; Delete; }

// Bottom signal path cutout in bottom ground
Box(8) = {-padW/2,-padW/2,-gndT, padW/2 + traceL+gap, padW, gndT};
BooleanDifference{ Volume{5}; }{ Volume{8}; Delete; }

// Top signal metal: pad on TSV + trace + gap + ground ext
Box(9)  = {-padW/2,-padW/2,h, padW, padW, traceT}; // pad
Box(10) = {padW/2,-traceW/2,h, traceL-padW/2, traceW, traceT}; // trace
Box(11) = {traceL,-traceW/2,h, gap, traceW, traceT}; // port 1 gap
Box(12) = {traceL+gap,-traceW/2,h, box-(traceL+gap), traceW, traceT}; // ground ext

// Bottom signal metal
Box(13) = {-padW/2,-padW/2,-gndT-traceT, padW, padW, traceT}; // pad
Box(14) = {padW/2,-traceW/2,-gndT-traceT, traceL-padW/2, traceW, traceT};
Box(15) = {traceL,-traceW/2,-gndT-traceT, gap, traceW, traceT}; // port 2 gap
Box(16) = {traceL+gap,-traceW/2,-gndT-traceT, box-(traceL+gap), traceW, traceT};
// Connection from TSV bottom to bottom pad
Box(17) = {-padW/2,-padW/2,-gndT, padW, padW, gndT}; // fills hole in ground z=-2..0

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,6,9,10,11,12,13,14,15,16,17}; Delete; };

air()={}; si()={}; ox()={}; cu()={}; gnd()={}; metal()={}; p1v()={}; p2v()={};
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  inSiZ = (zmin >= -0.1 && zmax <= h+0.1 && dz > 10);
  inBotGndZ = (zmax <= 0.1 && dz <= gndT+0.5);
  inTopGndZ = (zmin >= h-0.1 && zmax <= h+gndT+0.1 && dz <= gndT+0.5);
  inTopTrZ = (zmin >= h-0.1 && zmax <= h+traceT+0.1 && dz > gndT+0.5 || (zmin >= h-0.1 && zmax <= h+traceT+0.1 && dz > 0.5 && Abs(cy) <= padW/2+1));
  inBotTrZ = (zmin >= -gndT-traceT-0.1 && zmax <= -gndT+0.1 && dz > 0.5);
  isOuterXY = (dx > 20 || dy > 20);

  If (inBotGndZ)
    inPadArea = (Abs(cx) <= padW/2+1 && Abs(cy) <= padW/2+1 && dx <= padW+2);
    If (inPadArea) metal() += {v(i)};
    Else gnd() += {v(i)};
    EndIf
  ElseIf (inTopGndZ && isOuterXY)
    gnd() += {v(i)};
  ElseIf (inSiZ && !isOuterXY)
    If (dx <= tsv+0.2 && dy <= tsv+0.2) cu() += {v(i)};
    Else ox() += {v(i)};
    EndIf
  ElseIf (inSiZ && isOuterXY)
    si() += {v(i)};
  ElseIf (zmin >= h-0.1 && zmax <= h+traceT+0.1)
    isP1 = (Abs(cx - (traceL+gap/2)) < 2.0 && dx < gap+1.5 && Abs(cy) < traceW);
    If (isP1) p1v() += {v(i)};
    Else metal() += {v(i)};
    EndIf
  ElseIf (inBotTrZ)
    isP2 = (Abs(cx - (traceL+gap/2)) < 2.0 && dx < gap+1.5 && Abs(cy) < traceW);
    If (isP2) p2v() += {v(i)};
    Else metal() += {v(i)};
    EndIf
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
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.4;
For i In {0:#bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps || zmin < -airBot-gndT-traceT+eps || zmax > h+airTop+eps);
  isYFace = (dy < 0.5);
  inGapX1 = (xmin >= traceL-eps && xmax <= traceL+gap+eps);
  inGapZ1 = (zmin >= h-eps && zmax <= h+traceT+eps);
  inGapZ2 = (zmin >= -gndT-traceT-eps && zmax <= -gndT+eps);
  isP1 = (isYFace && inGapX1 && inGapZ1 && Abs(cy) <= traceW/2+eps);
  isP2 = (isYFace && inGapX1 && inGapZ2 && Abs(cy) <= traceW/2+eps);
  If (onOuter) absb() += {s};
  ElseIf (isP1) p1() += {s};
  ElseIf (isP2) p2() += {s};
  EndIf
EndFor
// PEC: boundaries of ground volumes only (not material interfaces, not ports)
For i In {0:#gnd()-1}
  gb() = Boundary{ Volume{gnd(i)}; };
  For j In {0:#gb()-1}
    gs = Abs(gb(j));
    bb() = BoundingBox Surface{gs};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
    dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
    onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps || zmin < -airBot-gndT-traceT+eps || zmax > h+airTop+eps);
    isYFace = (dy < 0.5);
    inGapX1 = (xmin >= traceL-eps && xmax <= traceL+gap+eps);
    inGapZ1 = (zmin >= h-eps && zmax <= h+traceT+eps);
    inGapZ2 = (zmin >= -gndT-traceT-eps && zmax <= -gndT+eps);
    isP1s = (isYFace && inGapX1 && inGapZ1 && Abs(cy) <= traceW/2+eps);
    isP2s = (isYFace && inGapX1 && inGapZ2 && Abs(cy) <= traceW/2+eps);
    If (!onOuter && !isP1s && !isP2s) pec() += {gs}; EndIf
  EndFor
EndFor
If (#p1() == 0 || #p2() == 0)
  allS() = Surface{:};
  For i In {0:#allS()-1}
    ss = allS(i);
    bb2() = BoundingBox Surface{ss};
    xmin2=bb2(0); ymin2=bb2(1); zmin2=bb2(2); xmax2=bb2(3); ymax2=bb2(4); zmax2=bb2(5);
    cx2=(xmin2+xmax2)/2; cy2=(ymin2+ymax2)/2; cz2=(zmin2+zmax2)/2;
    dx2=xmax2-xmin2; dy2=ymax2-ymin2; dz2=zmax2-zmin2;
    isY2 = (dy2 < 0.6);
    inGapX1b = (xmin2 >= traceL-0.5 && xmax2 <= traceL+gap+0.5);
    inGapZ1b = (zmin2 >= h-0.5 && zmax2 <= h+traceT+0.5);
    inGapZ2b = (zmin2 >= -gndT-traceT-0.5 && zmax2 <= -gndT+0.5);
    isP1b = (isY2 && inGapX1b && inGapZ1b && Abs(cy2) <= traceW/2+0.5);
    isP2b = (isY2 && inGapX1b && inGapZ2b && Abs(cy2) <= traceW/2+0.5);
    If (isP1b) p1() += {ss}; EndIf
    If (isP2b) p2() += {ss}; EndIf
  EndFor
EndIf
If (#absb() > 0) Physical Surface(10) = {absb()}; EndIf
If (#pec() > 0) Physical Surface(11) = {pec()}; EndIf
If (#p1() > 0) Physical Surface(12) = {p1()}; EndIf
If (#p2() > 0) Physical Surface(13) = {p2()}; EndIf
Mesh.CharacteristicLengthMin = 2.0; Mesh.CharacteristicLengthMax = 12.0; Mesh.Algorithm3D = 10;
