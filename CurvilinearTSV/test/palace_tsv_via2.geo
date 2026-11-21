SetFactory("OpenCASCADE");
// Units: microns (Palace L0=1e-6)
// Defensible 2-port TSV via: Cu TSV through SiO2/Si, bottom PEC ground,
// top PEC ground plane with signal pad. Lumped ports at top and bottom gaps.

tsv_w = 10.0; t_ox = 0.5; h_si = 50.0;
box = 60.0; air_top = 30.0; air_bot = 20.0; gnd_t = 2.0;
pad = 20.0; pad_t = 2.0; gap = 2.0;

// Domain
Box(1) = {-box,-box,-air_bot-gnd_t, 2*box,2*box, h_si+air_top+air_bot+gnd_t};
// Si
Box(2) = {-box,-box,0, 2*box,2*box,h_si};
// SiO2
Box(3) = {-tsv_w/2-t_ox,-tsv_w/2-t_ox,0, tsv_w+2*t_ox,tsv_w+2*t_ox,h_si};
// Cu TSV
Box(4) = {-tsv_w/2,-tsv_w/2,0, tsv_w,tsv_w,h_si};
// Bottom ground
Box(5) = {-box,-box,-gnd_t, 2*box,2*box,gnd_t};
// Top ground plane (covers top except signal pad region)
Box(6) = {-box,-box,h_si, 2*box,2*box,pad_t};
// Top signal pad (on top of TSV)
Box(7) = {-pad/2,-pad/2,h_si, pad,pad,pad_t};
// Bottom signal pad (under TSV)
Box(8) = {-pad/2,-pad/2,-gnd_t, pad,pad,gnd_t};
// Top port gap ring: cut gap around signal pad in top ground
// We'll model gap as a box ring by fragmenting; easier to make explicit gap volume
Box(9)  = { pad/2, -pad/2-gap, h_si, gap, pad+2*gap, pad_t}; // +x gap
Box(10) = {-pad/2-gap, -pad/2-gap, h_si, gap, pad+2*gap, pad_t}; // -x gap
Box(11) = {-pad/2,  pad/2, h_si, pad, gap, pad_t}; // +y gap
Box(12) = {-pad/2, -pad/2-gap, h_si, pad, gap, pad_t}; // -y gap
// Bottom port gap ring
Box(13) = { pad/2, -pad/2-gap, -gnd_t, gap, pad+2*gap, gnd_t};
Box(14) = {-pad/2-gap, -pad/2-gap, -gnd_t, gap, pad+2*gap, gnd_t};
Box(15) = {-pad/2,  pad/2, -gnd_t, pad, gap, gnd_t};
Box(16) = {-pad/2, -pad/2-gap, -gnd_t, pad, gap, gnd_t};

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:16}; Delete; };

// Classify
air()={}; si()={}; ox()={}; cu()={}; gnd()={}; metal()={}; portTop()={}; portBot()={};
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;

  inSiZ = (zmin >= -0.1 && zmax <= h_si+0.1 && dz > 10);
  inBotZ = (zmax <= 0.1);
  inTopZ = (zmin >= h_si-0.1 && zmax <= h_si+pad_t+0.1);
  isOuterXY = (dx > 25 || dy > 25);

  If (inBotZ && dz <= gnd_t+0.1)
    // bottom: ground, or port gap, or signal pad
    inPad = (Abs(cx) <= pad/2+0.1 && Abs(cy) <= pad/2+0.1 && dx <= pad+0.1 && dy <= pad+0.1);
    inGap = (!inPad && Abs(cx) <= pad/2+gap+0.1 && Abs(cy) <= pad/2+gap+0.1);
    If (inPad) metal() += {v(i)};
    ElseIf (inGap) portBot() += {v(i)};
    Else gnd() += {v(i)};
    EndIf
  ElseIf (inTopZ)
    inPad = (Abs(cx) <= pad/2+0.1 && Abs(cy) <= pad/2+0.1 && dx <= pad+0.1 && dy <= pad+0.1);
    inGap = (!inPad && Abs(cx) <= pad/2+gap+0.1 && Abs(cy) <= pad/2+gap+0.1);
    If (inPad) metal() += {v(i)};
    ElseIf (inGap) portTop() += {v(i)};
    Else gnd() += {v(i)};
    EndIf
  ElseIf (inSiZ && !isOuterXY)
    If (dx <= tsv_w+0.2 && dy <= tsv_w+0.2) cu() += {v(i)};
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
If (#portTop() > 0) Physical Volume(7) = {portTop()}; EndIf
If (#portBot() > 0) Physical Volume(8) = {portBot()}; EndIf

allv()={}; allv()+=air(); allv()+=si(); allv()+=ox(); allv()+=cu(); allv()+=gnd(); allv()+=metal(); allv()+=portTop(); allv()+=portBot();
bnd()={}; For i In {0:#allv()-1} bnd() += Boundary{ Volume{allv(i)}; }; EndFor

pec()={}; p1()={}; p2()={}; absb()={};
eps=0.3;
For i In {0:#bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps || zmin < -air_bot-gnd_t+eps || zmax > h_si+air_top+eps);
  // Port 1: +x face of top port gap at x=pad/2 (signal pad edge)
  isP1 = (Abs(cx - pad/2) < eps && dx < eps+0.1 && zmin >= h_si-eps && zmax <= h_si+pad_t+eps && Abs(cy) <= pad/2+eps);
  isP2 = (Abs(cx - pad/2) < eps && dx < eps+0.1 && zmax <= 0.1 && zmin >= -gnd_t-eps && Abs(cy) <= pad/2+eps);
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

Mesh.CharacteristicLengthMin = 1.5;
Mesh.CharacteristicLengthMax = 10.0;
Mesh.Algorithm3D = 10;
