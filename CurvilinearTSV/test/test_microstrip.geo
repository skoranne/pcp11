SetFactory("OpenCASCADE");
// um, L0=1e-6
box=80; hsub=10; air=40; gndT=2; traceW=4; traceT=2; traceL=60; gap=2;
Box(1) = {-box,-box,-gndT, 2*box,2*box, hsub+air+gndT}; // air
Box(2) = {-box,-box,0, 2*box,2*box,hsub}; // substrate
Box(3) = {-box,-box,-gndT, 2*box,2*box,gndT}; // ground
// trace: left gnd ext, gap1, signal, gap2, right gnd ext
Box(4) = {-traceL/2-gap-traceL/4,-traceW/2,hsub, traceL/4,traceW,traceT}; // left gnd
Box(5) = {-traceL/2-gap,-traceW/2,hsub, gap,traceW,traceT}; // port1 gap
Box(6) = {-traceL/2,-traceW/2,hsub, traceL,traceW,traceT}; // signal
Box(7) = {traceL/2,-traceW/2,hsub, gap,traceW,traceT}; // port2 gap
Box(8) = {traceL/2+gap,-traceW/2,hsub, traceL/4,traceW,traceT}; // right gnd
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,6,7,8}; Delete; };
airv()={}; subv()={}; gndv()={}; metal()={}; p1v()={}; p2v()={};
For i In {0:#v()-1}
  bb()=BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin; cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  If (zmax <= 0.1) gndv() += {v(i)};
  ElseIf (zmin >= -0.1 && zmax <= hsub+0.1 && dz > 1) subv() += {v(i)};
  ElseIf (zmin >= hsub-0.1)
    If (Abs(cx - (-traceL/2-gap/2)) < 1.5) p1v() += {v(i)};
    ElseIf (Abs(cx - (traceL/2+gap/2)) < 1.5) p2v() += {v(i)};
    Else metal() += {v(i)};
    EndIf
  Else airv() += {v(i)};
  EndIf
EndFor
If(#airv()>0) Physical Volume(1)={airv()}; EndIf
If(#subv()>0) Physical Volume(2)={subv()}; EndIf
If(#gndv()>0) Physical Volume(3)={gndv()}; EndIf
If(#metal()>0) Physical Volume(4)={metal()}; EndIf
If(#p1v()>0) Physical Volume(5)={p1v()}; EndIf
If(#p2v()>0) Physical Volume(6)={p2v()}; EndIf
allv()={}; allv()+=airv(); allv()+=subv(); allv()+=gndv(); allv()+=metal(); allv()+=p1v(); allv()+=p2v();
bnd()={}; For i In {0:#allv()-1} bnd()+=Boundary{Volume{allv(i)};}; EndFor
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.3;
For i In {0:#bnd()-1}
  s=Abs(bnd(i)); bb()=BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin<-box+eps || xmax>box-eps || ymin<-box+eps || ymax>box-eps || zmin<-gndT-eps || zmax>hsub+air-eps);
  isYFace = (dy < 0.5);
  inZ = (zmin >= hsub-eps && zmax <= hsub+traceT+eps);
  isP1 = (isYFace && Abs(cx-(-traceL/2-gap/2))<1.5 && inZ && Abs(cy)<=traceW/2+eps);
  isP2 = (isYFace && Abs(cx-(traceL/2+gap/2))<1.5 && inZ && Abs(cy)<=traceW/2+eps);
  If(onOuter) absb()+={s}; ElseIf(isP1) p1()+={s}; ElseIf(isP2) p2()+={s}; EndIf
EndFor
For i In {0:#gndv()-1}
  gb()=Boundary{Volume{gndv(i)};};
  For j In {0:#gb()-1}
    gs=Abs(gb(j)); bb()=BoundingBox Surface{gs};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    onOuter=(xmin<-box+eps||xmax>box-eps||ymin<-box+eps||ymax>box-eps||zmin<-gndT-eps||zmax>hsub+air-eps);
    If(!onOuter) pec()+={gs}; EndIf
  EndFor
EndFor
If(#absb()>0) Physical Surface(10)={absb()}; EndIf
If(#pec()>0) Physical Surface(11)={pec()}; EndIf
If(#p1()>0) Physical Surface(12)={p1()}; EndIf
If(#p2()>0) Physical Surface(13)={p2()}; EndIf
Mesh.CharacteristicLengthMin=1.5; Mesh.CharacteristicLengthMax=10; Mesh.Algorithm3D=10;
