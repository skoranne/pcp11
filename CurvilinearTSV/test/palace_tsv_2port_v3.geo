SetFactory("OpenCASCADE");
// 2-port TSV. Units: um, L0=1e-6. P1: +X edge top trace. P2: -X edge bottom trace.
tsv=10.0; tox=0.5; h=50.0; box=80.0; airTop=40.0; airBot=20.0; gndT=2.0;
traceW=6.0; traceT=2.0; padW=12.0;

Box(1) = {-box,-box,-airBot-gndT, 2*box,2*box, h+airTop+airBot+gndT}; // air
Box(2) = {-box,-box,0, 2*box,2*box,h}; // Si
Box(3) = {-tsv/2-tox,-tsv/2-tox,0, tsv+2*tox,tsv+2*tox,h}; // SiO2
Box(4) = {-tsv/2,-tsv/2,0, tsv,tsv,h}; // Cu TSV
Box(5) = {-box,-box,-gndT, 2*box,2*box,gndT}; // bottom ground
Box(6) = {-box,-box,h, 2*box,2*box,gndT};      // top ground
// Top cutout +X
Box(7) = {-padW/2,-padW/2,h, padW/2+box, padW, gndT};
BooleanDifference{ Volume{6}; }{ Volume{7}; Delete; }
// Bottom cutout -X
Box(8) = {-box,-padW/2,-gndT, box+padW/2, padW, gndT};
BooleanDifference{ Volume{5}; }{ Volume{8}; Delete; }
// Top signal: pad + trace to +X
Box(9)  = {-padW/2,-padW/2,h, padW, padW, traceT};
Box(10) = {padW/2,-traceW/2,h, box-padW/2, traceW, traceT};
// Bottom signal: pad + trace to -X
Box(11) = {-padW/2,-padW/2,-gndT-traceT, padW, padW, traceT};
Box(12) = {-box,-traceW/2,-gndT-traceT, box-padW/2, traceW, traceT};
// TSV-bottom connection
Box(13) = {-padW/2,-padW/2,-gndT, padW, padW, gndT};

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,6,9,10,11,12,13}; Delete; };

air()={}; si()={}; ox()={}; cu()={}; gnd()={}; metal()={};
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  inSiZ = (zmin >= -0.1) * (zmax <= h+0.1) * (dz > 10);
  inBotGndZ = (zmin >= -gndT-0.1) * (zmax <= 0.1) * (dz <= gndT+0.5);
  inTopGndZ = (zmin >= h-0.1) * (zmax <= h+gndT+0.1) * (dz <= gndT+0.5);
  inTopTrZ  = (zmin >= h-0.1) * (zmax <= h+traceT+0.1) * (dz > 0.5);
  inBotTrZ  = (zmin >= -gndT-traceT-0.1) * (zmax <= -gndT+0.1) * (dz > 0.5);
  isOuterXY = (dx > 20) + (dy > 20);
  inPadArea = (Abs(cx) <= padW/2+1) * (Abs(cy) <= padW/2+1) * (dx <= padW+2);

  If (inBotGndZ)
    If (inPadArea) metal() += {v(i)}; EndIf
    If (!inPadArea) gnd() += {v(i)}; EndIf
  EndIf
  If (inTopGndZ * isOuterXY)
    gnd() += {v(i)};
  EndIf
  If (inSiZ * (!isOuterXY))
    If (dx <= tsv+0.2 && dy <= tsv+0.2) cu() += {v(i)}; EndIf
    If (!(dx <= tsv+0.2 && dy <= tsv+0.2)) ox() += {v(i)}; EndIf
  EndIf
  If (inSiZ * isOuterXY)
    si() += {v(i)};
  EndIf
  If (inTopTrZ + inBotTrZ)
    metal() += {v(i)};
  EndIf
  If (!(inBotGndZ) * !(inTopGndZ * isOuterXY) * !(inSiZ) * !(inTopTrZ + inBotTrZ))
    air() += {v(i)};
  EndIf
EndFor

If(#air()>0) Physical Volume(1)={air()}; EndIf
If(#si()>0)  Physical Volume(2)={si()};  EndIf
If(#ox()>0)  Physical Volume(3)={ox()};  EndIf
If(#cu()>0)  Physical Volume(4)={cu()};  EndIf
If(#gnd()>0) Physical Volume(5)={gnd()}; EndIf
If(#metal()>0) Physical Volume(6)={metal()}; EndIf

allv()={}; allv()+=air(); allv()+=si(); allv()+=ox(); allv()+=cu(); allv()+=gnd(); allv()+=metal();
bnd()={}; For i In {0:#allv()-1} bnd()+=Boundary{Volume{allv(i)};}; EndFor
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.4;
For i In {0:#bnd()-1}
  s=Abs(bnd(i)); bb()=BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin < -box+eps) + (xmax > box-eps) + (ymin < -box+eps) + (ymax > box-eps) + (zmin < -airBot-gndT-traceT+eps) + (zmax > h+airTop+eps);
  onOuter = (onOuter > 0);
  isXFace = (dx < 0.5);
  inPortY = (Abs(cy) <= traceW/2+eps);
  inFullZ = (zmin <= h+traceT+5) * (zmax >= -gndT-traceT-5);
  isP1 = isXFace * (xmax > box-eps) * inPortY * inFullZ;
  isP2 = isXFace * (xmin < -box+eps) * inPortY * inFullZ;
  If (isP1) p1()+={s}; EndIf
  If (isP2) p2()+={s}; EndIf
  If (onOuter * (!isP1) * (!isP2)) absb()+={s}; EndIf
EndFor
For i In {0:#gnd()-1}
  gb()=Boundary{Volume{gnd(i)};};
  For j In {0:#gb()-1}
    gs=Abs(gb(j)); bb()=BoundingBox Surface{gs};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
    onOuter = (xmin < -box+eps) + (xmax > box-eps) + (ymin < -box+eps) + (ymax > box-eps) + (zmin < -airBot-gndT-traceT+eps) + (zmax > h+airTop+eps);
    onOuter = (onOuter > 0);
    isXFace=(dx<0.5); inPortY=(Abs(cy)<=traceW/2+eps);
    inFullZ=(zmin <= h+traceT+5)*(zmax >= -gndT-traceT-5);
    isP1s=isXFace*(xmax>box-eps)*inPortY*inFullZ;
    isP2s=isXFace*(xmin<-box+eps)*inPortY*inFullZ;
    If ((!onOuter)*(!isP1s)*(!isP2s)) pec()+={gs}; EndIf
  EndFor
EndFor
If(#absb()>0) Physical Surface(10)={absb()}; EndIf
If(#pec()>0) Physical Surface(11)={pec()}; EndIf
If(#p1()>0) Physical Surface(12)={p1()}; EndIf
If(#p2()>0) Physical Surface(13)={p2()}; EndIf
Mesh.CharacteristicLengthMin=2.0; Mesh.CharacteristicLengthMax=15.0; Mesh.Algorithm3D=1;
Field[1] = Box;
Field[1].VIn = 3.0; Field[1].VOut = 15;
Field[1].XMin = -box-1; Field[1].XMax = box+1;
Field[1].YMin = -traceW-2; Field[1].YMax = traceW+2;
Field[1].ZMin = -gndT-traceT-5; Field[1].ZMax = h+traceT+5;
Background Field = 1;
